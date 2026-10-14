"""
LAPSE: Linear Algebraic Probing Security Estimator (GF(2), GF(p), GF(2^m))

DSL quick reference
-------------------
Types:
  - `bit` or `gf2`                 → GF(2)
  - `gf(p)`                        → prime field, e.g. `gf(257)`
  - `gf2m(m, poly)`                → GF(2^m) with irreducible polynomial `poly` (integer), e.g.
                                     `gf2m(8, 0x11B)` → AES field (x^8 + x^4 + x^3 + x + 1)

Ops (linear-only):
  - `+` field addition (always allowed)
  - `-` field subtraction (same as `+` in char-2 fields)
  - `xor` only allowed when characteristic is 2 (alias of `+` there)
  - `*` scalar multiplication **only if one side is a numeric literal**

Examples:
  input gf(257) a[2]; const gf(257) c = 5; output gf(257) b[2];
  b[0] = a[0] + 3 * a[1] + c;             // linear over GF(257)

  input gf2m(8, 0x11B) x[2]; rand gf2m(8, 0x11B) r;
  y[0] = x[0] xor r;                       // ok (char-2)
  y[1] = x[1] + (0x02 * r);                // ok scalar mult in GF(2^8)

CLI:
  python -m lapse --demo_gf2 --t 1
  python -m lapse --demo_gfp --t 1
  python -m lapse --demo_gf2m --t 1

Notes:
- All declarations in one program must use the **same type** (field); mixed fields are rejected.
- `S` (the secret target) is never probeable.
- Random-probing δ(p) supports each of these finite fields.

Requires: lark  (pip install lark)
License: MIT
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Union
import itertools
import argparse
import sys
import random
import math

# =============== Finite fields ===============
class Field:
    def zero(self) -> int: raise NotImplementedError
    def one(self) -> int: raise NotImplementedError
    def add(self, a:int,b:int) -> int: raise NotImplementedError
    def sub(self, a:int,b:int) -> int: return self.add(a, self.neg(b))
    def neg(self, a:int) -> int: raise NotImplementedError
    def mul(self, a:int,b:int) -> int: raise NotImplementedError
    def inv(self, a:int) -> int: raise NotImplementedError
    def from_int(self, x:int) -> int: raise NotImplementedError
    def eq0(self, a:int) -> bool: return a == self.zero()
    @property
    def char(self) -> int: raise NotImplementedError

class PrimeField(Field):
    def __init__(self, p:int):
        if not isinstance(p, int) or p < 2 or p >= 2**64:
            raise ValueError("p must be a prime integer with 2 <= p < 2**64")
        # These Miller--Rabin bases are deterministic for unsigned 64-bit integers.
        if p not in (2, 3):
            if p % 2 == 0:
                raise ValueError("p must be prime")
            d, s = p - 1, 0
            while d % 2 == 0:
                d //= 2
                s += 1
            for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
                if a % p == 0:
                    continue
                x = pow(a, d, p)
                if x in (1, p - 1):
                    continue
                for _ in range(s - 1):
                    x = x * x % p
                    if x == p - 1:
                        break
                else:
                    raise ValueError("p must be prime")
        self.p = p
    def zero(self): return 0
    def one(self): return 1 % self.p
    @property
    def char(self): return self.p
    def from_int(self, x:int) -> int: return x % self.p
    def add(self,a,b): return (a+b) % self.p
    def neg(self,a): return (-a) % self.p
    def mul(self,a,b): return (a*b) % self.p
    def inv(self,a):
        if a % self.p == 0: raise ZeroDivisionError("inv(0)")
        # Fermat little theorem since p is prime
        return pow(a, self.p-2, self.p)

class GF2m(Field):
    def __init__(self, m:int, poly:int):
        if not isinstance(m, int) or m < 1:
            raise ValueError("m must be an integer >= 1")
        if not isinstance(poly, int) or poly < 0 or poly >> m != 1:
            # poly must be degree m (leading 1 at bit m)
            raise ValueError("poly must be degree m irreducible (leading 1 at bit m)")
        self.m = m
        self.poly = poly
        self.mask = (1<<m)-1
        # A degree-m polynomial is irreducible iff x^(2^m) = x modulo it
        # and it has no irreducible factor of degree at most floor(m/2).
        x = self._poly_remainder(2, poly)
        power = x
        for degree in range(1, m + 1):
            power = self.mul(power, power)
            if degree <= m // 2:
                a, b = power ^ x, poly
                while b:
                    a, b = b, self._poly_remainder(a, b)
                if a != 1:
                    raise ValueError("poly must be irreducible over GF(2)")
        if power != x:
            raise ValueError("poly must be irreducible over GF(2)")

    @staticmethod
    def _poly_remainder(a, b):
        while a and a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        return a
    @property
    def char(self): return 2
    def zero(self): return 0
    def one(self): return 1
    def from_int(self, x:int) -> int: return abs(x) & self.mask
    def add(self,a,b): return (a ^ b) & self.mask
    def neg(self,a): return a & self.mask  # char-2
    def mul(self,a,b):
        a &= self.mask; b &= self.mask
        res = 0
        x = a
        y = b
        for _ in range(self.m*2):
            if y & 1:
                res ^= x
            y >>= 1
            carry = x & (1<<(self.m-1))
            x = ((x<<1) & self.mask)
            if carry:
                x ^= (self.poly & self.mask)
            if y == 0: break
        return res & self.mask
    def inv(self,a):
        a &= self.mask
        if a == 0: raise ZeroDivisionError("inv(0)")
        exponent, result = (1 << self.m) - 2, 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.mul(a, a)
            exponent >>= 1
        return result

# =============== Parser ===============
try:
    from lark import Lark, Transformer, v_args, Token
except ImportError:
    print("This script requires 'lark' (pip install lark).", file=sys.stderr)
    raise

DSL_GRAMMAR = r"""
?start: program

program: CIRCUIT IDENT "(" [param_list] ")" ";" decl_or_stmt*
param_list: IDENT ("," IDENT)*

?decl_or_stmt: decl | stmt

// Terminals with priority (avoid IDENT eating keywords)
ROLE.2: "input" | "output" | "rand" | "temp" | "const"
BIT.2:  "bit" | "gf2"
CIRCUIT.2: "circuit"
GF.2:  "gf"
GF2M.2: "gf2m"
TO.2: "to"

// Declarations
role: ROLE

type: BIT                -> t_bit
    | GF "(" NUMBER ")"   -> t_gfp
    | GF2M "(" NUMBER "," NUMBER ")" -> t_gf2m

ident_arr: IDENT ("[" expr "]")?

decl: role type ident_arr ("=" expr)? ";"

?stmt: assign | for_stmt | call_stmt
call_stmt: call ";"
assign: lvalue "=" expr ";"
lvalue: IDENT ("[" expr "]")?

# allow optional "step EXPR" before ')'

for_stmt: "for" IDENT "in" "[" expr ((".."|TO) expr) ("step" expr)? ")" "{" decl_or_stmt* "}"



?expr: expr "+" term  -> addb
     | expr "-" term  -> subb
     | expr "xor" term -> xorb
     | term

?term: term "*" factor -> mulb
     | factor

?factor: NUMBER           -> number
       | IDENT ("[" expr "]")?   -> var
       | call
       | "(" expr ")"

call: IDENT "(" [arg_list] ")"
arg_list: expr ("," expr)*

%import common.CNAME -> IDENT
NUMBER: /[+-]?(?:0[xX][0-9a-fA-F]+|0[bB][01]+|[0-9]+)/
%import common.WS
%ignore WS
%ignore /\/\/[^\n]*/         // // comments
%ignore /\/\*.*?\*\//s      // /* block comments */
"""

# =============== AST ===============
@dataclass
class Program:
    name: str
    params: List[str]
    body: List[Union['Decl','Assign','For']]

@dataclass
class CallStmt:
    name: str
    args: List['Expr']

@dataclass
class Decl:
    role: str
    typeinfo: Tuple[str, Tuple[int, ...]]  # ('gf2',()), ('gfp',(p,)), ('gf2m',(m,poly))
    name: str
    size_expr: Optional['Expr'] = None
    init: Optional['Expr'] = None

@dataclass
class Assign:
    name: str
    index: Optional['Expr']
    expr: 'Expr'

@dataclass
class For:
    var: str
    start: 'Expr'
    end: 'Expr'
    step: Optional['Expr']  # None → step=1
    body: List[Union['Decl','Assign','For','CallStmt']]

class Expr: ...

@dataclass
class Number(Expr): value: int
@dataclass
class Var(Expr): name: str; index: Optional[Expr] = None
@dataclass
class Call(Expr): name: str; args: List[Expr]
@dataclass
class BinOp(Expr): op: str; left: Expr; right: Expr  # 'add','sub','xor','mul'

# =============== Transformer ===============
class ToAST(Transformer):
    @v_args(inline=False)
    def program(self, items):
        """
        Grammar: program: CIRCUIT IDENT "(" [param_list] ")" ";" decl_or_stmt*
        After child transforms:
        - Tokens (CIRCUIT, punctuation, etc.)
        - Token('IDENT', <name>)
        - param_list -> Python list[str] (optional)
        - decl_or_stmt* -> AST nodes
        """
        name = None
        params = []
        body = []

        for x in items:
            # take the first IDENT token as name
            if isinstance(x, Token) and x.type == 'IDENT' and name is None:
                name = x.value
                continue
            # The optional parameter list has already been transformed to strings.
            if isinstance(x, list) and not params:
                params = x
                continue
            # collect only AST nodes for the body
            if isinstance(x, (Decl, Assign, For, CallStmt)):
                body.append(x)

        if name is None:
            raise EvalError("Malformed program: missing circuit name")

        return Program(name=name, params=params, body=body)

    @v_args(inline=True)
    def param_list(self, *ids):
        return [t.value for t in ids]

    # type variants
    @v_args(inline=True)
    def t_bit(self, tok):
        return ('gf2', ())
    @v_args(inline=True)
    def t_gfp(self, _gf, p):
        return ('gfp', (_parse_integer(p.value),))
    @v_args(inline=True)
    def t_gf2m(self, _tag, m, poly):
        return ('gf2m', (_parse_integer(m.value), _parse_integer(poly.value)))

    @v_args(inline=True)
    def role(self, tok:Token): return tok

    @v_args(inline=True)
    def decl(self, role_tok, typeinfo, ident_arr, init=None):
        name, size_expr = ident_arr
        return Decl(role=role_tok.value, typeinfo=typeinfo, name=name, size_expr=size_expr, init=init)

    @v_args(inline=True)
    def ident_arr(self, ident, size=None):
        return ident.value, size

    @v_args(inline=True)
    def assign(self, lvalue, expr):
        name, idx = lvalue
        return Assign(name=name, index=idx, expr=expr)

    @v_args(inline=True)
    def lvalue(self, ident, idx=None):
        return ident.value, idx

    @v_args(inline=False)
    def for_stmt(self, items):
        """
        items may contain tokens '..' / 'to' and 'step'. Extract:
        IDENT, start_expr, end_expr, optional step_expr, and body nodes.
        """
        if not items or not isinstance(items[0], Token) or items[0].type != 'IDENT':
            raise EvalError("Malformed for_stmt: missing loop variable")

        var_tok = items[0]
        idx = 1

        def next_expr():
            nonlocal idx
            while idx < len(items) and isinstance(items[idx], Token):
                idx += 1
            if idx >= len(items):
                raise EvalError("Malformed for_stmt: missing expression")
            expr = items[idx]
            idx += 1
            return expr

        start = next_expr()
        while idx < len(items) and isinstance(items[idx], Token):
            idx += 1
        end = next_expr()

        while idx < len(items) and isinstance(items[idx], Token):
            idx += 1
        step = None
        if idx < len(items) and isinstance(items[idx], Expr):
            step = items[idx]
            idx += 1

        body = [node for node in items[idx:] if not isinstance(node, Token)]
        return For(var=var_tok.value, start=start, end=end, step=step, body=body)
    @v_args(inline=True)
    def call_stmt(self, call):
        assert isinstance(call, Call)
        return CallStmt(name=call.name, args=call.args)


    @v_args(inline=True)
    def number(self, n): return Number(_parse_integer(n.value))

    @v_args(inline=True)
    def var(self, ident, idx=None): return Var(ident.value, idx)

    @v_args(inline=True)
    def call(self, ident, args=None):
        if args is None: args=[]
        elif not isinstance(args, list): args=[args]
        return Call(ident.value, args)

    @v_args(inline=True)
    def arg_list(self, *args): return list(args)

    @v_args(inline=True)
    def addb(self, a, b): return BinOp('add', a, b)
    @v_args(inline=True)
    def subb(self, a, b): return BinOp('sub', a, b)
    @v_args(inline=True)
    def xorb(self, a, b): return BinOp('xor', a, b)
    @v_args(inline=True)
    def mulb(self, a, b): return BinOp('mul', a, b)

    def IDENT(self, tok:Token): return tok

# =============== Utilities ===============
class EvalError(Exception): pass

def _parse_integer(value: str) -> int:
    unsigned = value.lstrip('+-').lower()
    return int(value, 0 if unsigned.startswith(('0x', '0b')) else 10)

def eval_int(expr: Expr, params: Dict[str,int]) -> int:
    if isinstance(expr, Number): return expr.value
    if isinstance(expr, Var):
        if expr.index is not None: raise EvalError("Index expr cannot be array access")
        if expr.name not in params: raise EvalError(f"Unknown integer symbol {expr.name}")
        return int(params[expr.name])
    if isinstance(expr, BinOp):
        a = eval_int(expr.left, params)
        b = eval_int(expr.right, params)
        if expr.op == 'add': return a + b
        if expr.op == 'sub': return a - b
        if expr.op == 'xor': return a ^ b
        if expr.op == 'mul': return a * b
    if isinstance(expr, Call):
        if expr.name == 'pow2' and len(expr.args) == 1:
            exponent = eval_int(expr.args[0], params)
            if exponent < 0:
                raise EvalError("pow2 exponent must be non-negative")
            return 1 << exponent
        raise EvalError(f"Unsupported integer call {expr.name} in index")

    raise EvalError(f"Unsupported integer expr: {expr}")
    
# =============== IR/SSA ===============
@dataclass
class Symbol:
    role: str
    typeinfo: Tuple[str, Tuple[int,...]]
    size: Optional[int]

@dataclass
class SSAValue:
    name: str
    index: Optional[int]
    version: int
    role: str
    def ssa(self):
        return f"{self.name}{f'[{self.index}]' if self.index is not None else ''}#{self.version}"

@dataclass
class IR:
    symbols: Dict[str, Symbol]
    versions: Dict[Tuple[str,Optional[int]], int]
    field: Field
    # source columns & vectors
    sources: List[str]                 # list of base source SSA names (column order)
    vectors: Dict[str, List[int]]      # ssa -> vector over sources (field elements)
    created: List[Tuple[str,str]]      # (ssa, role)
    secret_ssa: Optional[str]

# =============== Builder ===============
class Builder:
    def __init__(self, params: Dict[str,int]):
        self.params = params.copy()
        self.symbols: Dict[str, Symbol] = {}
        self.versions: Dict[Tuple[str,Optional[int]], int] = {}
        self._initialized = set()
        self.field: Optional[Field] = None
        self._field_typeinfo = None
        self.sources: List[str] = []
        self.vectors: Dict[str, List[int]] = {}
        self.created: List[Tuple[str,str]] = []
        self.secret_ssa: Optional[str] = None
        self.rand_counter: int = 0  # fresh IDs for randbit() calls


    # ----- field selection -----
    def _ensure_field(self, typeinfo: Tuple[str, Tuple[int,...]]):
        if self.field is not None:
            if typeinfo != self._field_typeinfo:
                raise EvalError("All declarations must use the same field type")
            return
        tag, params = typeinfo
        if tag == 'gf2':
            self.field = GF2m(1, 0b11)  # GF(2) as m=1 with x+1 (any irreducible of deg 1)
        elif tag == 'gfp':
            (p,) = params
            self.field = PrimeField(p)
        elif tag == 'gf2m':
            m, poly = params
            self.field = GF2m(m, poly)
        else:
            raise EvalError(f"Unknown field tag {tag}")
        self._field_typeinfo = typeinfo

            # ----- symbols -----
    def declare(self, d: Decl):
        # 0) duplicate check up-front
        if d.name in self.symbols:
            raise EvalError(f"Duplicate declaration {d.name}")
        if d.name.startswith('__') or d.name == 'CONST_ONE':
            raise EvalError(f"Symbol name {d.name} is reserved for internal use")
        if d.name == 'S' and (d.role != 'const' or d.size_expr is not None):
            raise EvalError("The secret target S must be declared as a scalar const")
        if d.size_expr is not None and d.init is not None:
            raise EvalError("Array initializers are not supported")
        if d.role == 'input' and d.init is not None:
            raise EvalError("Inputs are independent sources and cannot have initializers")

        # 1) ensure field is chosen
        self._ensure_field(d.typeinfo)

        # 2) resolve size once
        size = None
        if d.size_expr is not None:
            size = eval_int(d.size_expr, self.params)
            if size < 0:
                raise EvalError("Array size must be non-negative")

        # 3) register the symbol BEFORE materializing anything
        self.symbols[d.name] = Symbol(role=d.role, typeinfo=d.typeinfo, size=size)

        # 4) materialize inputs eagerly (sources); others stay lazy
        if d.role == 'input':
            if size is None:
                self._materialize_input_source(d.name, None)
            else:
                for i in range(size):
                    self._materialize_input_source(d.name, i)

    # helper: add one input element as a fresh source column, widening existing rows
    def _materialize_input_source(self, name: str, idx: Optional[int]):
        ssa = self._next(name, idx, 'input')

        # register new source column
        src_name = ssa.ssa()
        self._as_source(src_name)  # must append to self.sources

        # widen all existing row vectors to new width with zeros
        new_C = len(self.sources)
        zero = self.field.zero()
        for k, row in self.vectors.items():
            if len(row) < new_C:
                row.extend([zero] * (new_C - len(row)))

        # create one-hot vector for this input source
        vec = [zero] * new_C
        vec[new_C - 1] = self.field.one()
        self.vectors[src_name] = vec
    
    # ----- SSA helpers -----
    def _next(self, name:str, index:Optional[int], role:str) -> SSAValue:
        key = (name, index)
        v = self.versions.get(key, -1) + 1
        self.versions[key] = v
        self._initialized.add(key)
        ssa = SSAValue(name, index, v, role)
        self.created.append((ssa.ssa(), role))
        return ssa

    def _as_source(self, ssa_name:str):
        if ssa_name not in self.sources:
            self.sources.append(ssa_name)
            # extend existing vectors with 0 in new column
            for vec in self.vectors.values():
                vec.append(self.field.zero())

    def _validate_access(self, name: str, index: Optional[int]) -> Symbol:
        sym = self.symbols.get(name)
        if sym is None:
            raise EvalError(f"Use of undeclared symbol {name}")
        if sym.size is None and index is not None:
            raise EvalError(f"Scalar {name} cannot be indexed")
        if sym.size is not None:
            if index is None:
                raise EvalError(f"Array {name} requires an index")
            if not 0 <= index < sym.size:
                raise EvalError(f"Index {index} outside {name}[0:{sym.size}]")
        return sym

    def _validate_range(self, name: str, start: int, end: int) -> Symbol:
        sym = self.symbols.get(name)
        if sym is None or sym.size is None:
            raise EvalError(f"Expected a declared array: {name}")
        if not 0 <= start <= end <= sym.size:
            raise EvalError(f"Invalid range [{start}, {end}) for {name}[{sym.size}]")
        return sym

    def _read_base(self, name:str, index:Optional[int]) -> str:
        sym = self._validate_access(name, index)
        key = (name, index)
        if key not in self._initialized:
            if sym.role not in ("input", "rand"):
                raise EvalError(f"Reading {name} before definition")
            ssa = self._next(name, index, sym.role)
            self._as_source(ssa.ssa())
            vec = [self.field.zero() for _ in self.sources]
            vec[self.sources.index(ssa.ssa())] = self.field.one()
            self.vectors[ssa.ssa()] = vec
        return f"{name}{f'[{index}]' if index is not None else ''}#{self.versions[key]}"

    # ----- expression linearization -----
    def _vec_const(self, k:int) -> List[int]:
        # Public affine offsets do not affect recoverability. Only coefficients
        # of unknown input shares and randomness belong in the source matrix.
        return [self.field.zero() for _ in self.sources]

    def _vec_add(self, a: List[int], b: List[int]) -> List[int]:
        f = self.field
        # Pad to the current number of sources to keep vectors aligned
        L = max(len(a), len(b), len(self.sources))
        if len(a) < L:
            a = a + [f.zero()] * (L - len(a))
        if len(b) < L:
            b = b + [f.zero()] * (L - len(b))
        return [f.add(x, y) for x, y in zip(a, b)]


    def _vec_scale(self, k:int, v:List[int]) -> List[int]:
        f = self.field
        kk = f.from_int(k)
        return [f.mul(kk, x) for x in v]

    def _vec_sub(self, a: List[int], b: List[int]) -> List[int]:
        return self._vec_add(a, [self.field.neg(x) for x in b])

    def _vec_of(self, expr: Expr) -> List[int]:
        f = self.field

        if isinstance(expr, Call):
            if expr.name in ('rand', 'randbit'):
                if expr.args:
                    raise EvalError(f"{expr.name}() takes no arguments")
                # default the field if not set yet
                if self.field is None:
                    self.field = GF2m(1, 0b11)  # GF(2)

                # 1) allocate a new source column (no SSA row!)
                col_name = f"__randcall_{self.rand_counter}"
                self.rand_counter += 1
                self._as_source(col_name)

                # 2) widen all existing row vectors to the new width
                new_C = len(self.sources)
                zero = self.field.zero()
                for k, row in self.vectors.items():
                    if len(row) < new_C:
                        row.extend([zero] * (new_C - len(row)))

                # 3) return the one-hot vector on the new source column
                vec = [zero] * new_C
                vec[new_C - 1] = self.field.one()
                return vec

        if isinstance(expr, Number):
            return self._vec_const(expr.value)
      
        elif isinstance(expr, Var):
            if expr.index is None:
                ssa = self._read_base(expr.name, None)
            else:
                idx = eval_int(expr.index, self.params)
                ssa = self._read_base(expr.name, idx)
            return list(self.vectors[ssa])
        
        if isinstance(expr, BinOp):
            if expr.op in ('add','sub','xor'):
                va = self._vec_of(expr.left)
                vb = self._vec_of(expr.right)
                if expr.op == 'xor' and f.char != 2:
                    raise EvalError("'xor' allowed only in characteristic 2 fields")
                if expr.op == 'sub' and f.char == 2:
                    # same as add
                    return self._vec_add(va, vb)
                if expr.op == 'sub':
                    return self._vec_sub(va, vb)
                return self._vec_add(va, vb)
            if expr.op == 'mul':
                # allow only scalar * vector (one side literal)
                la = isinstance(expr.left, Number)
                lb = isinstance(expr.right, Number)
                if la and lb:
                    return self._vec_const(0)
                if la and not lb:
                    vb = self._vec_of(expr.right)
                    return self._vec_scale(expr.left.value, vb)
                if lb and not la:
                    va = self._vec_of(expr.left)
                    return self._vec_scale(expr.right.value, va)
                raise EvalError("'*' allowed only for scalar*expr (one side must be a literal)")
        if isinstance(expr, Call):
            if expr.name == 'xor_reduce' or expr.name == 'sum_reduce':
                if len(expr.args) != 3 or not isinstance(expr.args[0], Var) or expr.args[0].index is not None:
                    raise EvalError("xor_reduce/sum_reduce requires (array_name, start, end)")
                arr = expr.args[0].name
                s = eval_int(expr.args[1], self.params)
                e = eval_int(expr.args[2], self.params)
                self._validate_range(arr, s, e)
                if expr.name == 'xor_reduce' and f.char != 2:
                    raise EvalError("xor_reduce allowed only in characteristic 2 fields")
                acc = [self.field.zero() for _ in self.sources]  # may change length later; keep consistent
                # to be robust against new sources created during loop, recompute size if needed
                first = True
                for i in range(s, e):
                    ssa = self._read_base(arr, i)
                    vec = self.vectors[ssa]
                    if first:
                        acc = [self.field.zero() for _ in self.sources]
                        first = False
                    acc = self._vec_add(acc, vec)
                return acc
            raise EvalError(f"Unknown function call {expr.name}")
        raise EvalError(f"Unsupported expression: {expr}")

    # ----- build passes -----
    def build(self, prog: Program):
        self._execute_block(prog.body)

    def _execute_block(self, nodes):
        # Register declarations first so input source order is deterministic.
        # Initializers are executable statements and run exactly once, in order.
        for node in nodes:
            if isinstance(node, Decl):
                self.declare(node)
        for node in nodes:
            if isinstance(node, Decl):
                if node.init is not None:
                    self._assign(Assign(node.name, None, node.init))
            elif isinstance(node, Assign):
                self._assign(node)
            elif isinstance(node, For):
                self._for(node)
            elif isinstance(node, CallStmt):
                self._exec_callstmt(node)
            else:
                raise EvalError(f"Unknown statement: {type(node)}")
    
    def _for(self, loop: For):
        s = eval_int(loop.start, self.params)
        e = eval_int(loop.end,   self.params)
        k = eval_int(loop.step,  self.params) if loop.step is not None else 1
        if k == 0:
            raise EvalError("for-loop step cannot be 0")
        previous = self.params.get(loop.var)
        had_previous = loop.var in self.params
        local_names = {node.name for node in loop.body if isinstance(node, Decl)}
        if local_names.intersection(self.symbols):
            raise EvalError("Loop-local declarations cannot shadow existing symbols")
        if any(isinstance(node, Decl) and (node.role in ('input', 'output') or node.name == 'S')
               for node in loop.body):
            raise EvalError("Inputs, outputs, and the target S must be declared at top level")
        try:
            for i in range(s, e, k):
                self.params[loop.var] = i
                try:
                    self._execute_block(loop.body)
                finally:
                    for name in local_names:
                        self.symbols.pop(name, None)
                    self._initialized.difference_update(
                        key for key in list(self._initialized) if key[0] in local_names
                    )
        finally:
            if had_previous:
                self.params[loop.var] = previous
            else:
                self.params.pop(loop.var, None)

    
    def _exec_callstmt(self, cs: CallStmt):
        if cs.name == 'refresh_range':
            # args: (IDENT array, start, end)
            if not (len(cs.args) == 3 and isinstance(cs.args[0], Var) and cs.args[0].index is None):
                raise EvalError("refresh_range(arr, start, end) requires a bare array name and two integer bounds")
            arr = cs.args[0].name
            L   = eval_int(cs.args[1], self.params)
            R   = eval_int(cs.args[2], self.params)
            sym = self._validate_range(arr, L, R)
            if sym.role == 'const':
                raise EvalError("Cannot refresh a const array")
            if R <= L: return
            hub = R - 1
            for i in range(L, R-1):
                # fresh r
                r_vec = self._vec_of(Call('randbit', []))
                r_ssa = self._next("__r", None, "rand")
                self.vectors[r_ssa.ssa()] = list(r_vec)
                # arr[i] += r; arr[hub] -= r, retaining every created SSA wire.
                ssa_i = self._read_base(arr, i); vec_i = self.vectors[ssa_i]
                ssa_h = self._read_base(arr, hub); vec_h = self.vectors[ssa_h]
                out_i = self._next(arr, i, sym.role)
                out_h = self._next(arr, hub, sym.role)
                self.vectors[out_i.ssa()] = self._vec_add(vec_i, r_vec)
                self.vectors[out_h.ssa()] = self._vec_sub(vec_h, r_vec)
            return
        raise EvalError(f"Unknown statement call {cs.name}")

    
    def _assign(self, node):
        # ---- extract RHS (supports expr|rhs|value) ----
        rhs = getattr(node, 'expr', None) or getattr(node, 'rhs', None) or getattr(node, 'value', None)
        if rhs is None:
            raise EvalError(f"Assign missing RHS; attrs={dir(node)}")

        # compute RHS vector BEFORE touching LHS version
        rhs_vec = self._vec_of(rhs)

        # ---- extract LHS in either style ----
        name = None
        idx_raw = None

        # Style A: lvalue/Var node
        lv = getattr(node, 'lvalue', None) or getattr(node, 'lhs', None) or getattr(node, 'target', None)
        if lv is not None:
            if not isinstance(lv, Var):
                raise EvalError("Assignment lvalue must be a Var")
            name = lv.name
            idx_raw = lv.index

        # Style B: flattened fields name/index
        if name is None and hasattr(node, 'name'):
            name = node.name
            idx_raw = getattr(node, 'index', None)

        if name is None:
            raise EvalError(f"Assign missing LHS; attrs={dir(node)}")

        # ---- evaluate optional index safely ----
        def _eval_opt_index(x):
            if x is None:
                return None
            # already int?
            if isinstance(x, int):
                return x
            # expression from our AST?
            if isinstance(x, Expr):
                return eval_int(x, self.params)
            # try generic int conversion (e.g., Token)
            try:
                return int(x)
            except Exception:
                raise EvalError(f"Unsupported index type for {name}: {type(x)}")

        idx = _eval_opt_index(idx_raw)

        # ---- look up symbol / role, bump version, store vector ----
        sym = self._validate_access(name, idx)
        if sym.role == 'const' and (name, idx) in self._initialized:
            raise EvalError(f"Cannot reassign const {name}")
        if sym.role == 'const' and name != 'S' and any(not self.field.eq0(x) for x in rhs_vec):
            raise EvalError(f"const {name} must be public; use temp for a derived wire")
        out = self._next(name, idx, sym.role)       # e.g., a[0]#1
        self.vectors[out.ssa()] = rhs_vec

        # ---- mark the secret when writing const S (scalar) ----
        if name == 'S' and sym.role == 'const' and idx is None:
            self.secret_ssa = out.ssa()

# =============== Linear algebra over general field ===============
def _latest_output_ssa_per_index(ir: IR) -> Dict[Tuple[str,int], str]:
        """
        For each output array element name[i], pick the SSA with the largest version.
        Returns {(name, i) -> ssa_name}.
        """
        latest: Dict[Tuple[str,int], Tuple[int,str]] = {}
        for ssa, role in ir.created:
            if role != 'output':
                continue
            # parse name[index]#version
            name, rest = ssa.split('#', 1)
            ver = int(rest)
            if '[' in name and name.endswith(']'):
                base = name[:name.rfind('[')]
                idx  = int(name[name.rfind('[')+1:-1])
            else:
                base, idx = name, None
            key = (base, idx)
            if key not in latest or ver > latest[key][0]:
                latest[key] = (ver, ssa)
        return {k:v for k,(_,v) in latest.items()}

def _input_ssa_per_index(ir: IR) -> Dict[Tuple[str,int], str]:
        """First (base) SSA per input index (usually version 0)."""
        base: Dict[Tuple[str,int], str] = {}
        for ssa, role in ir.created:
            if role != 'input':
                continue
            name, rest = ssa.split('#', 1)
            if '[' in name and name.endswith(']'):
                base_name = name[:name.rfind('[')]
                idx       = int(name[name.rfind('[')+1:-1])
            else:
                base_name, idx = name, None
            key = (base_name, idx)
            # keep earliest version encountered
            if key not in base:
                base[key] = ssa
        return base

def rank_over_field(rows: List[List[int]], F: Field) -> int:
        """Row rank via Gaussian elimination over a field, without mutating rows."""
        if not rows:
            return 0
        # make a copy; pad rows to equal width
        C = max((len(r) for r in rows), default=0)
        A = [r[:] + [F.zero()] * (C - len(r)) for r in rows]
        R = len(A)
        r = c = 0
        rank = 0
        while r < R and c < C:
            piv = None
            for i in range(r, R):
                if not F.eq0(A[i][c]):
                    piv = i; break
            if piv is None:
                c += 1; continue
            if piv != r:
                A[r], A[piv] = A[piv], A[r]
            inv = F.inv(A[r][c])
            for j in range(c, C):
                A[r][j] = F.mul(A[r][j], inv)
            for i in range(R):
                if i == r or F.eq0(A[i][c]):
                    continue
                factor = A[i][c]
                for j in range(c, C):
                    A[i][j] = F.sub(A[i][j], F.mul(factor, A[r][j]))
            r += 1; c += 1; rank += 1
        return rank

def candidates_inputs_outputs(ir: IR) -> Tuple[List[str], List[str]]:
        """Return (inputs, final_outputs) as SSA name lists."""
        inputs  = list(_input_ssa_per_index(ir).values())
        outputs = list(_latest_output_ssa_per_index(ir).values())
        return inputs, outputs

# ========= GF(2) fast path: packed elimination =========
try:
    import numpy as _np
except Exception:
    _np = None

def _pad_rows_to_width_gf2(rows, target_len):
    """Pad each row with zeros (ints) to target_len, coerce entries to {0,1}."""
    out = []
    for r in rows:
        rr = list(r) + [0] * (target_len - len(r))
        # coerce to 0/1 in case elements are field objects/ints
        out.append([int(x) & 1 for x in rr])
    return out

def _pack_uint64_gf2(rows_01):
    """Pack a 0/1 matrix into uint64 words: shape (R, W), W = ceil(C/64)."""
    R = len(rows_01)
    C = len(rows_01[0]) if R else 0
    W = (C + 63) // 64
    M = _np.zeros((R, W), dtype=_np.uint64)
    for i, row in enumerate(rows_01):
        for c, v in enumerate(row):
            if v & 1:
                w = c >> 6
                b = c & 63
                M[i, w] |= (_np.uint64(1) << _np.uint64(b))
    return M, C, W

def _rank_gf2_packed_numpy(M, C):
    """
    Row rank over GF(2) using packed uint64 matrix M (shape R x W).
    Column-wise elimination; complexity ~ O(R * C * W/64).
    """
    R, W = M.shape
    rank = 0
    r = 0
    for c in range(C):
        w = c >> 6
        mask = _np.uint64(1) << _np.uint64(c & 63)
        # find pivot at/under r with bit set in column c
        piv = -1
        # vectorized: which rows have the bit?
        colbits = (M[r:, w] & mask) != 0
        if colbits.any():
            piv = r + int(_np.argmax(colbits))
        if piv < 0:
            continue
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        # eliminate other rows having bit c
        # rows with bit set (excluding pivot row)
        rows_with_bit = (M[:, w] & mask) != 0
        rows_with_bit[r] = False
        idxs = _np.nonzero(rows_with_bit)[0]
        if idxs.size:
            M[idxs, :] ^= M[_np.newaxis, r, :]
        r += 1
        rank += 1
        if r == R:
            break
    return rank

def _in_rowspan_gf2_packed(rows, target):
    """
    Return True iff target is in the row-span of rows over GF(2).
    Uses NumPy packed elimination if available; else pure Python.
    """
    # Normalize widths and coerce to 0/1
    C = max(len(target), max((len(r) for r in rows), default=0))
    rows01 = _pad_rows_to_width_gf2(rows, C)
    targ01 = _pad_rows_to_width_gf2([target], C)[0]

    if _np is not None:
        M, C, W = _pack_uint64_gf2(rows01 + [targ01])
        # compute rank before/after adding target row
        rank_rows = _rank_gf2_packed_numpy(M[:-1].copy(), C)
        rank_aug  = _rank_gf2_packed_numpy(M.copy(), C)
        return rank_aug == rank_rows
    else:
        # ---- Pure Python bitset fallback (no NumPy) ----
        # pack each row as Python int (LSB = column 0)
        def pack_int(bits):
            v = 0
            for i, b in enumerate(bits):
                if b & 1:
                    v |= (1 << i)
            return v
        mat = [pack_int(r) for r in rows01]
        tgt = pack_int(targ01)
        # Gaussian elimination on integers
        rank = 0
        r = 0
        for c in range(C):
            mask = 1 << c
            piv = -1
            for i in range(r, len(mat)):
                if mat[i] & mask:
                    piv = i; break
            if piv < 0:
                continue
            if piv != r:
                mat[r], mat[piv] = mat[piv], mat[r]
            # eliminate from all other rows
            for i in range(len(mat)):
                if i != r and (mat[i] & mask):
                    mat[i] ^= mat[r]
            # track rank
            if mat[r]:
                rank += 1
                r += 1
                if r == len(mat):
                    break
        # After elimination, target is in span iff it reduces to 0 when augmented;
        # simulate by re-running with target included and comparing ranks:
        mat2 = [pack_int(r) for r in rows01] + [tgt]
        rank2 = 0
        r2 = 0
        for c in range(C):
            mask = 1 << c
            piv = -1
            for i in range(r2, len(mat2)):
                if mat2[i] & mask:
                    piv = i; break
            if piv < 0:
                continue
            if piv != r2:
                mat2[r2], mat2[piv] = mat2[piv], mat2[r2]
            for i in range(len(mat2)):
                if i != r2 and (mat2[i] & mask):
                    mat2[i] ^= mat2[r2]
            if mat2[r2]:
                rank2 += 1
                r2 += 1
                if r2 == len(mat2):
                    break
        return rank2 == rank


def _make_gf2_span_event(rows, target):
    """Prepack GF(2) rows once for repeated subset recovery queries.

    The returned function takes row indices. Each query constructs a triangular
    basis of the selected rows and reduces the target incrementally. Python's
    arbitrary-width integers avoid both repeated NumPy allocation and a second
    elimination of the augmented matrix.
    """
    def pack(row):
        value = 0
        for index, coefficient in enumerate(row):
            if int(coefficient) & 1:
                value |= 1 << index
        return value

    packed = [pack(row) for row in rows]
    packed_target = pack(target)

    def event(indices):
        remaining = packed_target
        if remaining == 0:
            return True
        basis = {}
        for index in indices:
            value = packed[index]
            while value:
                pivot = value.bit_length() - 1
                previous = basis.get(pivot)
                if previous is None:
                    basis[pivot] = value
                    break
                value ^= previous
            while remaining:
                previous = basis.get(remaining.bit_length() - 1)
                if previous is None:
                    break
                remaining ^= previous
            if remaining == 0:
                return True
        return False

    return event

def in_span(rows: List[List[int]], target: List[int], F: Field) -> bool:
    """
    Return True iff target is in the span of rows over field F.
    Fast GF(2) path via packed elimination; otherwise generic Gaussian elimination.
    """
    # ---- GF(2) fast path ----
    if (isinstance(F, PrimeField) and F.p == 2) or (isinstance(F, GF2m) and F.m == 1):
        # Extension fields also have characteristic two, but their elements
        # are not binary coefficients and require the generic field backend.
        return _in_rowspan_gf2_packed(rows, target)

    # ---- Generic field elimination ----
    if not rows:
        return all(F.eq0(x) for x in target)

    m_target = len(target)
    max_row = max((len(r) for r in rows), default=0)
    C = max(m_target, max_row)

    # pad each row and the target to width C
    A = [(r[:] + [F.zero()] * (C - len(r))) for r in rows]
    A.append(target[:] + [F.zero()] * (C - len(target)))
    R = len(A)

    r = c = 0
    while r < R and c < C:
        piv = None
        for i in range(r, R - 1):  # exclude target row from pivoting
            if not F.eq0(A[i][c]):
                piv = i; break
        if piv is None:
            c += 1
            continue
        if piv != r:
            A[r], A[piv] = A[piv], A[r]
        inv = F.inv(A[r][c])
        for j in range(c, C):
            A[r][j] = F.mul(A[r][j], inv)
        for i in range(R):
            if i == r or F.eq0(A[i][c]):
                continue
            factor = A[i][c]
            for j in range(c, C):
                A[i][j] = F.sub(A[i][j], F.mul(factor, A[r][j]))
        r += 1
        c += 1

    # target row zero ⇒ in span
    return all(F.eq0(x) for x in A[R - 1])

# =============== Checks ===============
@dataclass
class CheckResult:
    ok: bool
    reason: str
    counterexample: Optional[List[str]] = None


def filter_leakables(created: List[Tuple[str,str]], secret_ssa: Optional[str], allow_roles: Optional[set]) -> List[str]:
    """Return SSA leakage sites; every version of the analysis target is excluded."""
    target_name = secret_ssa.split('#', 1)[0] if secret_ssa is not None else 'S'
    if allow_roles is not None:
        unknown = set(allow_roles) - {'input', 'output', 'rand', 'temp', 'const'}
        if unknown:
            raise ValueError(f"Unknown leakage roles: {', '.join(sorted(unknown))}")
    return [ssa for ssa, role in created
            if ssa.split('#', 1)[0] != target_name
            and (role != 'const' if allow_roles is None else role in allow_roles)]


def t_probing(ir: IR, t: int, allow_roles: Optional[set]) -> CheckResult:
    _validate_count(t, 't', minimum=0)
    if ir.secret_ssa is None:
        return CheckResult(False, "No secret S defined.")
    S = ir.vectors[ir.secret_ssa]
    if all(ir.field.eq0(value) for value in S):
        return CheckResult(False, "Target S is public; the empty probe set reveals it.", [])
    leaks = filter_leakables(ir.created, ir.secret_ssa, allow_roles)
    rows = { s: ir.vectors[s] for s in leaks }
    # search minimal k ≤ t
    for k in range(1, min(t, len(leaks)) + 1):
        for comb in itertools.combinations(leaks, k):
            if in_span([rows[s] for s in comb], S, ir.field):
                return CheckResult(False, f"Found probe set of size {k} revealing S.", list(comb))
    return CheckResult(True, f"PASS: No probe set of size ≤ {t} reveals S.")

# Random probing E_X(p) / delta(p)

def _validate_probability(value: float, name: str, *, interior: bool = False) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite probability") from exc
    valid = math.isfinite(value) and (0 < value < 1 if interior else 0 <= value <= 1)
    if not valid:
        interval = "(0, 1)" if interior else "[0, 1]"
        raise ValueError(f"{name} must be finite and in {interval}")
    return value


def _validate_count(value: int, name: str, *, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _recovery_problem(ir: IR, allow_roles: Optional[set]):
    if ir.secret_ssa is None or ir.secret_ssa not in ir.vectors:
        raise ValueError("No target S defined; declare and assign const S before estimation")
    leaks = filter_leakables(ir.created, ir.secret_ssa, allow_roles)
    rows = [ir.vectors[s] for s in leaks]
    target = ir.vectors[ir.secret_ssa]
    known = all(ir.field.eq0(x) for x in target)

    if ((isinstance(ir.field, PrimeField) and ir.field.p == 2)
            or (isinstance(ir.field, GF2m) and ir.field.m == 1)):
        return leaks, _make_gf2_span_event(rows, target), known

    def event(indices):
        return known or in_span([rows[i] for i in indices], target, ir.field)

    return leaks, event, known


@dataclass(frozen=True)
class RecoveryEstimate:
    """Recovery probability and sampling diagnostics.

    IS intervals use a normal approximation and can be unreliable when a few
    weighted observations dominate. A zero-hit IS run has no reported interval.
    The raw IS estimate is deliberately not clipped to [0, 1]: clipping biases it.
    ``pilot_samples`` counts observations discarded before the main estimate.
    """
    estimate: float
    wire_count: int
    method: str
    samples: int
    hits: int
    variance: Optional[float]
    standard_error: Optional[float]
    confidence_interval: Optional[Tuple[float, float]]
    confidence_level: float = 0.95
    confidence_method: str = "exact"
    effective_sample_size: Optional[float] = None
    proposal_rate: Optional[float] = None
    seed: Optional[int] = None
    pilot_samples: int = 0
    warnings: Tuple[str, ...] = ()


class _LogSumExp:
    """Accumulate nonnegative weighted observations without forming large weights."""
    __slots__ = ("max_log", "sum_scaled")

    def __init__(self):
        self.max_log = -math.inf
        self.sum_scaled = 0.0

    def add(self, logw: float):
        if logw == -math.inf:
            return
        if not math.isfinite(logw):
            raise ArithmeticError("Non-finite importance weight")
        if self.max_log == -math.inf:
            self.max_log, self.sum_scaled = logw, 1.0
        elif logw > self.max_log:
            self.sum_scaled = self.sum_scaled * math.exp(self.max_log - logw) + 1.0
            self.max_log = logw
        else:
            self.sum_scaled += math.exp(logw - self.max_log)

    def log_value(self) -> float:
        return self.max_log + math.log(self.sum_scaled) if self.sum_scaled else -math.inf

    def value(self) -> float:
        return _exp_or_inf(self.log_value())


def _exp_or_inf(log_value: float) -> float:
    try:
        return math.exp(log_value)
    except OverflowError:
        return math.inf


def _wilson_interval(hits: int, samples: int) -> Tuple[float, float]:
    z = 1.959963984540054
    mean = hits / samples
    denominator = 1.0 + z * z / samples
    center = (mean + z * z / (2 * samples)) / denominator
    half_width = z * math.sqrt(mean * (1 - mean) / samples + z * z / (4 * samples * samples)) / denominator
    return max(0.0, center - half_width), min(1.0, center + half_width)


def _sample_event(event, count: int, rate: float, samples: int, rng):
    rnd = rng.random
    for _ in range(samples):
        subset = [i for i in range(count) if rnd() < rate]
        yield bool(event(subset)), len(subset)


def _tune_proposal(event, count: int, p: float, budget: int, rng):
    """Pilot bisection aiming for a 5--30% event rate; discard every pilot draw."""
    lo, hi = p, 1.0
    q = max(p, 0.25)
    used = 0
    batch = max(1, budget // 5)
    while used < budget:
        size = min(batch, budget - used)
        hits = sum(hit for hit, _ in _sample_event(event, count, q, size, rng))
        used += size
        rate = hits / size
        if 0.05 <= rate <= 0.30:
            break
        if rate < 0.05:
            lo = q
        else:
            hi = q
        q = (lo + hi) / 2.0
    # Guard rounding at the endpoints; target p itself is never altered.
    return min(q, math.nextafter(1.0, 0.0)), used


def estimate_event(
    event,
    wire_count: int,
    p: float,
    *,
    method: str = "auto",
    samples: int = 10000,
    max_enum: int = 20,
    q: Optional[float] = None,
    seed: Optional[int] = None,
    rng=None,
    pilot_samples: int = 1000,
) -> RecoveryEstimate:
    """Estimate a Boolean event under independent Bernoulli wire leakage.

    ``event(indices)`` receives the indices of leaked sites, from 0 through
    ``wire_count - 1``. Conditional analyses must remove conditioned-absent sites
    before constructing the callback and wire count. The callback must be
    deterministic; monotonicity is not required. Diagnostics and pilot/main
    sampling follow :func:`estimate_recovery`; exact enumeration uses all sites.
    Automatic proposal tuning is a heuristic for nonmonotone event probabilities.
    """
    _validate_count(wire_count, 'wire_count', minimum=0)
    return _estimate_probability(
        event, wire_count, p, method=method, samples=samples, max_enum=max_enum,
        q=q, seed=seed, rng=rng, pilot_samples=pilot_samples)


def estimate_recovery(
    ir: IR,
    p: float,
    allow_roles: Optional[set] = None,
    *,
    method: str = "auto",
    samples: int = 10000,
    max_enum: int = 20,
    q: Optional[float] = None,
    seed: Optional[int] = None,
    rng=None,
    pilot_samples: int = 1000,
) -> RecoveryEstimate:
    """Estimate E_X with exact enumeration, MC or Bernoulli importance sampling.

    ``auto`` enumerates up to ``max_enum`` wires; otherwise an independent pilot
    selects MC or IS, followed by fresh main observations. When IS has no explicit
    ``q``, an additional pilot tunes it. ``samples`` is the main-run budget only.
    A supplied ``rng`` must expose ``random()``; use ``seed`` or ``rng``, not both.
    Endpoint probabilities and a public/zero target are evaluated deterministically.
    Identical vectors and zero rows are integrated out during exact enumeration.
    """
    from dataclasses import replace

    leaks, event, known = _recovery_problem(ir, allow_roles)
    # Duplicate-vector sites matter only through whether any representative leaks.
    groups = {}
    for index, ssa in enumerate(leaks):
        row = list(ir.vectors[ssa])
        while row and ir.field.eq0(row[-1]):
            row.pop()
        if not row:
            continue
        key = tuple(row)
        if key in groups:
            groups[key][1] += 1
        else:
            groups[key] = [index, 1]
    result = _estimate_probability(
        event, len(leaks), p, method=method, samples=samples, max_enum=max_enum,
        q=q, seed=seed, rng=rng, pilot_samples=pilot_samples,
        known=known, exact_groups=list(groups.values()))
    if result.method == 'exact':
        return result
    input_sites = set(_input_ssa_per_index(ir).values())
    input_indices = [i for i, ssa in enumerate(leaks) if ssa in input_sites]
    if input_indices and event(input_indices):
        lower_bound = float(p) ** len(input_indices)
        if result.estimate < lower_bound:
            warning = (
                f"Estimate is below the guaranteed input-subset recovery bound {lower_bound:.6g}; "
                "sampling has not resolved the recovery probability reliably."
            )
            result = replace(result, warnings=result.warnings + (warning,))
    return result


def _estimate_probability(
    event,
    count: int,
    p: float,
    *,
    method: str = "auto",
    samples: int = 10000,
    max_enum: int = 20,
    q: Optional[float] = None,
    seed: Optional[int] = None,
    rng=None,
    pilot_samples: int = 1000,
    known: bool = False,
    exact_groups=None,
) -> RecoveryEstimate:
    """Shared estimator engine; optional exact groups are recovery-specific."""
    p = _validate_probability(p, "p")
    _validate_count(samples, "samples")
    _validate_count(max_enum, "max_enum", minimum=0)
    _validate_count(pilot_samples, "pilot_samples")
    if q is not None:
        q = _validate_probability(q, "q", interior=True)
    if method not in {"auto", "exact", "mc", "is"}:
        raise ValueError("method must be one of auto, exact, mc, is")
    if seed is not None and rng is not None:
        raise ValueError("Specify either seed or rng, not both")
    rng = random.Random(seed) if rng is None else rng

    def exact_result(value):
        return RecoveryEstimate(value, count, "exact", 0, 0, 0.0, 0.0,
                                (value, value), seed=seed)

    if known:
        return exact_result(1.0)
    if count == 0 or p == 0.0:
        return exact_result(float(event([])))
    if p == 1.0:
        return exact_result(float(event(range(count))))
    if method == "exact" or (method == "auto" and count <= max_enum):
        if count > max_enum:
            raise ValueError(f"Exact enumeration needs 2^{count} patterns; max_enum={max_enum}")
        groups = exact_groups if exact_groups is not None else [(i, 1) for i in range(count)]
        representatives = [index for index, _ in groups]
        miss = [math.exp(multiplicity * math.log1p(-p)) for _, multiplicity in groups]
        present = [-math.expm1(multiplicity * math.log1p(-p)) for _, multiplicity in groups]
        group_count = len(representatives)
        terms = (
            math.prod(present[i] if (mask >> i) & 1 else miss[i]
                      for i in range(group_count))
            for mask in range(1 << group_count)
            if event([representatives[i] for i in range(group_count) if (mask >> i) & 1])
        )
        return exact_result(math.fsum(terms))

    pilot_used = 0
    if method == "auto":
        pilot_hits = sum(hit for hit, _ in _sample_event(event, count, p, pilot_samples, rng))
        pilot_used += pilot_samples
        method = "is" if pilot_hits < 10 else "mc"
    if method == "is" and q is None:
        q, used = _tune_proposal(event, count, p, pilot_samples, rng)
        pilot_used += used
    proposal = p if method == "mc" else q
    warnings = []
    if method == "mc":
        hits = sum(hit for hit, _ in _sample_event(event, count, p, samples, rng))
        mean = hits / samples
        variance = mean * (1 - mean) / (samples - 1) if samples > 1 else None
        se = math.sqrt(variance) if variance is not None else None
        interval = _wilson_interval(hits, samples)
        ess = float(hits)
        interval_method = "Wilson 95%"
        if hits == 0:
            warnings.append("No event occurrences observed; zero is not evidence of zero probability.")
    else:
        log_ratio_hit = math.log(p) - math.log(proposal)
        log_ratio_miss = math.log1p(-p) - math.log1p(-proposal)
        weights, squared = _LogSumExp(), _LogSumExp()
        hits = 0
        for hit, size in _sample_event(event, count, proposal, samples, rng):
            if hit:
                hits += 1
                log_weight = size * log_ratio_hit + (count - size) * log_ratio_miss
                weights.add(log_weight)
                squared.add(2.0 * log_weight)
        log_sum, log_sum2 = weights.log_value(), squared.log_value()
        mean = _exp_or_inf(log_sum - math.log(samples))
        ess = min(float(hits), _exp_or_inf(2 * log_sum - log_sum2)) if hits else 0.0
        if samples > 1:
            if hits:
                # Compute the variance and SE separately: squaring a tiny SE can underflow.
                log_fraction = min(0.0, 2 * log_sum - math.log(samples) - log_sum2)
                remainder = -math.expm1(log_fraction)
                log_var = (log_sum2 - math.log(samples) - math.log(samples - 1)
                           + math.log(remainder)) if remainder > 0 else -math.inf
                variance = _exp_or_inf(log_var)
                se = _exp_or_inf(log_var / 2.0)
            else:
                variance, se = 0.0, 0.0
        else:
            variance, se = None, None
        if hits and se is not None and math.isfinite(mean) and math.isfinite(se):
            z = 1.959963984540054
            interval = (max(0.0, min(1.0, mean - z * se)),
                        max(0.0, min(1.0, mean + z * se)))
            interval_method = "normal approximation 95% (importance sampling)"
        else:
            interval, interval_method = None, "unavailable"
        if hits == 0:
            warnings.append("No event occurrences under the proposal; IS uncertainty is unresolved. Increase samples or change q.")
        elif hits < 10 or ess < 10:
            warnings.append("Few effective event observations; the IS interval may be unreliable.")
        if hits and mean == 0.0:
            warnings.append("The estimate underflows binary64; use a less extreme parameter regime.")
        if not math.isfinite(mean) or (se is not None and not math.isfinite(se)):
            warnings.append("Importance weights exceed binary64 range; change q closer to p.")
    if samples < 2:
        warnings.append("At least two main samples are required to estimate variance.")
    return RecoveryEstimate(mean, count, method, samples, hits, variance, se, interval,
                            confidence_method=interval_method, effective_sample_size=ess,
                            proposal_rate=proposal, seed=seed, pilot_samples=pilot_used,
                            warnings=tuple(warnings))


# =============== Parsing and circuit construction ===============

def parse_program(text:str) -> Program:
    parser = Lark(DSL_GRAMMAR, parser='lalr', start='start')
    tree = parser.parse(text)
    return ToAST().transform(tree)

def build_ir(prog: Program, params: Dict[str,int]) -> IR:
    b = Builder(params.copy())
    if any(not isinstance(value, int) for value in params.values()):
        raise EvalError("Circuit parameters must be integers")
    b.build(prog)

    if b.field is None:
        b.field = GF2m(1, 0b11)

    return IR(symbols=b.symbols, versions=b.versions, field=b.field,
              sources=b.sources, vectors=b.vectors, created=b.created,
              secret_ssa=b.secret_ssa)

# =============== Demos ===============
DEMO_GF2 = """
// GF(2) two-share mask
circuit two_share_mask(n);
input bit a[2];
rand  bit r;
output bit b[2];

b[0] = a[0] xor r;
b[1] = a[1] xor r;

const bit S;
S = xor_reduce(a, 0, 2);
"""

DEMO_GFP = """
// GF(257) linear gadget
circuit gfp_demo(n);
input gf(257) a[2];
output gf(257) b[2];
const gf(257) c = 5;

b[0] = a[0] + 3 * a[1] + c;
b[1] = a[1] - a[0] + 2;   // linear over GF(257)

const gf(257) S;
S = a[0] + a[1];
"""

DEMO_GF2M = """
// GF(2^8) with AES polynomial 0x11B
circuit gf2m_demo(n);
input gf2m(8, 0x11B) x[2];
rand  gf2m(8, 0x11B) r;
output gf2m(8, 0x11B) y[2];

y[0] = x[0] xor r;
y[1] = x[1] + 0x02 * r;   // scalar mul

const gf2m(8, 0x11B) S;
S = x[0] + x[1];
"""

# =============== CLI ===============

def main(argv=None):
    """Analyze a circuit; use --json to retain statistics and provenance."""
    import json
    from pathlib import Path
    from lark.exceptions import LarkError
    from lapse.cli_utils import parse_params, parse_roles, provenance, result_dict, clean_json

    ap = argparse.ArgumentParser(description="LAPSE: Linear Algebraic Probing Security Estimator")
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument('--file', help='Path to .circ file')
    source.add_argument('--demo_gf2', action='store_true')
    source.add_argument('--demo_gfp', action='store_true')
    source.add_argument('--demo_gf2m', action='store_true')
    ap.add_argument('--n', type=int, help='Positive value for parameter n')
    ap.add_argument('--params', default='', help='Comma-separated integer parameters, e.g. rounds=3')
    ap.add_argument('--t', type=int, help='Check recovery of S by at most t probes')
    ap.add_argument('--allow-roles', help='Comma-separated probeable roles (default: exclude const)')
    ap.add_argument('--delta-p', type=float, help='Recovery probability E_X at probing rate p')
    ap.add_argument('--beta-p', type=float, help='Conditional boundary exposure at probing rate p < 1')
    ap.add_argument('--boundary-share', help='Boundary share for beta, e.g. a[0] or b[0] (or exact SSA name)')
    ap.add_argument('--beta-policy', choices=['separated', 'literal'], default='separated',
                    help='Explicit boundary reduction convention; see docs/beta-definition.md')
    ap.add_argument('--beta-q', type=float, help='IS proposal rate for beta; implies method is for beta')
    ap.add_argument('--method', choices=['auto', 'exact', 'mc', 'is'], default='auto')
    ap.add_argument('--max-enum', type=int, default=20, help='Maximum wire count for exact enumeration (default: 20)')
    ap.add_argument('--mc-samples', type=int, default=200000, help='Number of main MC/IS samples')
    ap.add_argument('--pilot-samples', type=int, default=1000, help='Independent pilot budget for automatic selection')
    ap.add_argument('--mc-q', type=float, help='IS proposal rate; implies --method is')
    ap.add_argument('--seed', type=int, default=0, help='Sampling seed (default: 0)')
    ap.add_argument('--json', action='store_true', help='Print machine-readable result and provenance')
    args = ap.parse_args(argv)
    try:
        if args.n is not None and args.n < 1:
            raise ValueError('n must be positive')
        if args.t is not None and args.t < 0:
            raise ValueError('t must be nonnegative')
        if args.mc_samples < 1 or args.pilot_samples < 1 or args.max_enum < 0:
            raise ValueError('sample counts must be positive and max-enum nonnegative')
        if args.mc_q is not None and args.method not in ('auto', 'is'):
            raise ValueError('--mc-q requires --method auto or is')
        if args.mc_q is not None and args.delta_p is None:
            raise ValueError('--mc-q requires --delta-p')
        if args.beta_p is not None and args.boundary_share is None:
            raise ValueError('--beta-p requires --boundary-share, e.g. a[0]')
        if args.beta_q is not None and (args.beta_p is None or args.method not in ('auto', 'is')):
            raise ValueError('--beta-q requires --beta-p and method auto or is')
        params = parse_params(args.params)
        if args.n is not None:
            if 'n' in params:
                raise ValueError('provide n only through --n or --params')
            params['n'] = args.n
        roles = parse_roles(args.allow_roles)
        text = (Path(args.file).read_text(encoding='utf-8') if args.file else
                DEMO_GFP if args.demo_gfp else DEMO_GF2M if args.demo_gf2m else DEMO_GF2)
        ir = build_ir(parse_program(text), params)
        output = {
            'field_characteristic': ir.field.char,
            'sources': len(ir.sources),
            'probeable_wires': len(filter_leakables(ir.created, ir.secret_ssa, roles)),
            'provenance': provenance(text, params, vars(args)),
        }
        if args.t is not None:
            from dataclasses import asdict
            output['t_probing'] = asdict(t_probing(ir, args.t, roles))
        if args.delta_p is not None:
            output['recovery'] = result_dict(estimate_recovery(
                ir, args.delta_p, roles, method='is' if args.mc_q is not None else args.method,
                samples=args.mc_samples, max_enum=args.max_enum, q=args.mc_q,
                seed=args.seed, pilot_samples=args.pilot_samples))
        if args.beta_p is not None:
            from lapse.boundary import estimate_beta
            output['beta'] = result_dict(estimate_beta(
                ir, args.beta_p, args.boundary_share, policy=args.beta_policy, allow_roles=roles,
                method='is' if args.beta_q is not None else args.method, samples=args.mc_samples,
                max_enum=args.max_enum, q=args.beta_q, seed=args.seed, pilot_samples=args.pilot_samples))
        if args.json:
            print(json.dumps(clean_json(output), indent=2, allow_nan=False))
        else:
            print(f"LAPSE: characteristic {ir.field.char}, {len(ir.sources)} sources, {output['probeable_wires']} probeable wires")
            if 't_probing' in output:
                check = output['t_probing']
                print(check['reason'])
                if check['counterexample']:
                    print('Counterexample: ' + ', '.join(check['counterexample']))
            if 'recovery' in output:
                result = output['recovery']
                print(f"E_X(p={args.delta_p:g}) = {result['estimate']:.12g} ({result['method']}, seed={args.seed})")
                print(f"Standard error: {result['standard_error']}; 95% interval: {result['confidence_interval']}")
                for warning in result['warnings']:
                    print(f"WARNING: {warning}")
            if 'beta' in output:
                result = output['beta']
                print(f"beta(p={args.beta_p:g}, U={result['boundary_share']}, policy={result['policy']}) = {result['estimate']:.12g}")
                print(f"Effective marginal rate p' = {result['p_effective']:.12g}")
                print(f"Standard error: {result['standard_error']}; 95% interval: {result['confidence_interval']}")
                for warning in result['warnings']:
                    print(f"WARNING: {warning}")
    except (ValueError, OSError, EvalError, LarkError) as exc:
        ap.error(str(exc))
