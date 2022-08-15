# SpatialPy is a Python 3 package for simulation of
# spatial deterministic/stochastic reaction-diffusion-advection problems
# Copyright (C) 2019 - 2022 SpatialPy developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
# published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU GENERAL PUBLIC LICENSE Version 3 for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import ast

class BuildExpression:
    """
    Accepts an expression string to validate and convert.
    Allows for pre-flight syntax and namespace validations,
    as well as converting between Python and C++ expressions.

    Object for managing context to validate Python expressions.
    Expressions can be passed and validated, which are validated for syntax, namespace, and other conditions.
    The later provided statements are expected to be valid Python expressions.

    :param blacklist: List of operators which are not allowed in the following expressions.
    Note that this will be "forwarded" to all following expressions.
    Ideally, one should define the "universal" blacklist in the constructor,
    using the `BuildExpression#with_blacklist` method for more granular validations.
    :type blacklist: list[str]

    :param namespace: Dictionary mapping allowed bare identifiers to their sanitized equivalents.
    Any bare identifiers not listed as a namespace key will trigger a failed validation.
    :type  namespace: dict[str, any]

    :param sanitize: Whether or not to substitute namespace names during conversion.
    Any valid names found as namespace keys will automatically be converted to the
    corresponding namespace values in the namespace dict when getexpr_* methods are called.
    :type sanitize: bool
    """
    def __init__(self, blacklist=None, namespace=None, sanitize=False):
        if blacklist is None:
            blacklist = dict({})
        elif not isinstance(blacklist, dict):
            blacklist = {ast_op: op for op, ast_op in BuildExpression.map_operator(blacklist)}
        if namespace is None:
            namespace = {}
        self.blacklist = blacklist
        self.namespace = namespace
        self.sanitize = sanitize

    class ValidationVisitor(ast.NodeTransformer):
        """
        A subclass of ast.NodeTransformer used to sanitize spatialpy expresions.

        :param namespace: SpactialPy namespace.
        :type namespace: dict

        :param blacklist: A list of blacklist operators.
        :type blacklist: list

        :param sanitize: Whether or not to sanitize nodes.
        :type sanitize: bool
        """
        def __init__(self, namespace=None, blacklist=None, sanitize=False):
            self.namespace = dict({}) if namespace is None else namespace
            self.blacklist = dict({}) if blacklist is None else blacklist
            self.sanitize = sanitize
            self.invalid_names = []
            self.invalid_operators = []

        def _check_blacklist(self, operator):
            """
            Helper function to check literal expression operators against blacklisted operators

            :param operator: operator to check
            :type operator: ast.BinOp.Op | ast.UnaryOp.Op
            """
            operator = type(operator)
            if operator in self.blacklist:
                self.invalid_operators.append(str(self.blacklist.get(operator)))

        def visit_Name(self, node): # pylint: disable=invalid-name
            """
            Visitor function for ast.Name nodes.

            :param node: A node of type 'Name' from an AST.
            :type node: ast.Name

            :returns: The original node if sanitize is False, otherwise a sanitized node.
            :rtype: ast.Name
            """
            if node.id not in self.namespace:
                self.invalid_names.append(node.id)
            elif self.sanitize:
                node.id = self.namespace.get(node.id)
            self.generic_visit(node)
            return node

        def visit_Call(self, node): # pylint: disable=invalid-name
            """
            Visitor function for ast.Call nodes.

            :param node: A node of type 'Call' from an AST.
            :type node: ast.Call

            :returns: The original node if sanitize is False, otherwise a sanitized node.
            :rtype: ast.Call
            """
            if node.func.id not in self.namespace:
                self.invalid_names.append(node.func.id)
            elif self.sanitize:
                node.func.id = self.namespace.get(node.func.id)
            self.generic_visit(node)
            return node

        def visit_BinOp(self, node): # pylint: disable=invalid-name
            """
            Visitor function for ast.BinOp nodes.

            :param node: A node of type 'BinOp' from an AST.
            :type node: ast.BinOp

            :returns: The original node.
            :rtype: ast.BinOp
            """
            self._check_blacklist(node.op)
            self.generic_visit(node)
            return node

        def visit_UnaryOp(self, node): # pylint: disable=invalid-name
            """
            Visitor function for ast.UnaryOp nodes.

            :param node: A node of type 'UnaryOp' from an AST.
            :type node: ast.UnaryOp

            :returns: The original node.
            :rtype: ast.UnaryOp
            """
            self._check_blacklist(node.op)
            self.generic_visit(node)
            return node

        def visit_BoolOp(self, node): # pylint: disable=invalid-name
            """
            Visitor function for ast.BoolOp nodes.

            :param node: A node of type 'BoolOp' from an AST.
            :type node: ast.BoolOp

            :returns: The original node.
            :rtype: ast.BinBoolOpOp
            """
            self._check_blacklist(node.op)
            self.generic_visit(node)
            return node

        def visit_Compare(self, node): # pylint: disable=invalid-name
            """
            Visitor function for ast.Compare nodes.

            :param node: A node of type 'Compare' from an AST.
            :type node: ast.Compare

            :returns: The original node.
            :rtype: ast.Compare
            """
            for operator in node.ops:
                self._check_blacklist(operator)
            self.generic_visit(node)
            return node

        def visit_Assign(self, node: "ast.Assign"): # pylint: disable=invalid-name
            """
            Visitor function for ast.Assign nodes.

            :param node: A node of type 'Assign' from an AST.
            :type node: ast.Assign

            :returns: The original node.
            :rtype: ast.Assign
            """
            self._check_blacklist(ast.Assign())
            self.generic_visit(node)
            return node

    operator_map = {
        # Basic math operators
        "+": ast.Add, "-": ast.Sub, "*": ast.Mult, "/": ast.Div,
        "**": ast.Pow, "//": ast.FloorDiv, "%": ast.Mod, "@": ast.MatMult,
        # Variable operators
        "=": ast.Assign, ":=": ast.Assign,
        # Boolean operators
        "and": ast.And, "or": ast.Or,
        # Bitwise operators (^ gets substituted for ** in Python and pow() in C++)
        "^": ast.BitXor, "<<": ast.LShift, ">>": ast.RShift, "|": ast.BitOr, "&": ast.BitAnd,
        # Comparison operators
        "==": ast.Eq, "!=": ast.NotEq, "!": ast.Not, ">": ast.Gt,
        ">=": ast.GtE, "<": ast.Lt, "<=": ast.LtE,
    }

    def __get_expr(self, converter):
        validator = BuildExpression.ValidationVisitor(self.namespace, self.blacklist, self.sanitize)
        validator.visit(converter.tree)

        if validator.invalid_operators:
            return None

        if validator.invalid_names:
            return None

        return converter.get_str()

    @classmethod
    def map_operator(cls, operator):
        """
        Map operator strings with built expressions.

        :param operator: Operator to be mapped
        :type operator: str | list[str]
        """
        if isinstance(operator, list):
            for oper in operator:
                yield from BuildExpression.map_operator(oper)
        else:
            # Base case: operator is a single string.
            if operator in BuildExpression.operator_map:
                yield operator, BuildExpression.operator_map.get(operator)
            elif operator in BuildExpression.operator_map.values():
                # Yield the operator directly if there is no need to map it.
                yield operator

    def with_blacklist(self, blacklist=None):
        """
        Create a new duplicate of the current expression, with a different operator blacklist.
        Overrides operator handling behavior when converting or validating the expression.

        :param blacklist: List of operators which are not allowed.
        :type blacklist: list[str]

        :returns: New expression containing the given blacklist.
        The returned expression is a *copy* of the current expression.
        :rtype: BuildExpression
        """
        if blacklist is None:
            blacklist = self.blacklist.copy()
        return BuildExpression(blacklist=blacklist, namespace=self.namespace)

    def with_namespace(self, namespace=None):
        """
        Create a new duplicate of the current expression, with a different namespace.
        Any identifiers present in the expression which are not listed in the namespace
        will cause the expression to be flagged as an invalid namespace during validation.

        :param namespace: A dictionary containing the namespace mappings for the expression.
        The keys of the dict are expected to be the "only" valid identifiers.
        The values of the namespace are what the keys map to during sanitization, if used.
        :type namespace: dict[str, str]

        :returns: New expression containing the given namespace.
        The returned expression is a *copy* of the current expression.
        :rtype: BuildExpression
        """
        if namespace is None:
            namespace = self.namespace.copy()
        return BuildExpression(blacklist=self.blacklist.copy(), namespace=namespace)

    def validate(self, statement):
        """
        Using the information provided so far, ensure that the given Python expression is valid.
        The Python expression is parsed, raising a SyntaxError if it is an invalid Python expression.
        The expression is then checked against the given properties, such as namespace and operator blacklist.
        Additionally, the expression is rejected if it is not a single rvalue expression.

        :param statement: A python expression.
        :type statement: str

        :returns: Result object containing the lists of invalid names and invalid operators.
        :rtype: ExpressionResults

        :raises SyntaxError: The statement is not a valid Python expression.
        """
        expr = ast.parse(statement)
        validator = BuildExpression.ValidationVisitor(self.namespace, self.blacklist, self.sanitize)
        validator.visit(expr)

        return ExpressionResults(invalid_names=validator.invalid_names, invalid_operators=validator.invalid_operators)

    def getexpr_python(self, statement):
        """
        Converts the expression object into a Python expression string.

        :param statement: C++ expression to be converted.
        :type statement: str

        :returns: Python expression string, if valid. Returns None if validation fails.
        :rtype: str | None

        :raises SyntaxError: If the C++ expression could not be converted to a Python expression.
        """
        expr = ast.parse(statement)
        return self.__get_expr(PythonConverter(expr))

    def getexpr_cpp(self, statement):
        """
        Converts the expression object into a C++ expression string.
        Raises a SyntaxError if conversion to a C++ string is impossible.

        :param statement: Python expression to be conveter.
        :type statement: str

        :returns: C++ expression string, if valid. Returns None if validation fails.
        :rtype: str | None

        :raises SyntaxError: If the Python expression could not be converted to a C++ expression.
        """
        statement = ExpressionConverter.convert_str(statement)
        expr = ast.parse(statement)
        return self.__get_expr(CppConverter(expr))


class ExpressionResults:
    """
    Container struct for returning the results of expression validation.
    Any expression items which indicate an invalid expression are listed on an ExpressionResults instance.
    Empty lists indicate that the expression is valid.

    Container struct for returning the results of expression validation.

    :param invalid_names: List of expression identifiers which were not valid in the given namespace.
    :type invalid_names: list[str]

    :param invalid_operators: List of blacklisted operators which were present in the expression.
    :type invalid_operators: list[str]

    :param is_valid: Override value for the `is_valid` property.
    If not set, then the validity of the expression is inferred by the `invalid_*` lists provided.
    :type is_valid: bool
    """
    def __init__(self, invalid_names=None, invalid_operators=None, is_valid=True):
        self.invalid_names = invalid_names
        self.invalid_operators = invalid_operators
        self.is_valid = is_valid and (not invalid_names and not invalid_operators)


class ExpressionConverter(ast.NodeVisitor):
    """
    A subclass of ast.NodeVisitor used to convert spatialpy expresions.

    :param tree: An abstract syntax tree.
    :type tree: ast.AST
    """
    def __init__(self, tree):
        self.tree = tree
        self.expression = []

    def _get_str(self, expr):
        self.visit(expr)
        return "".join(self.expression)

    @classmethod
    def convert_str(cls, expression):
        """
        Convert '^' to python pow operator.

        :param expression: BuildExpression to be converted.
        :type expression: str
        """
        return expression.replace("^", "**")

    def parse_operator(self, operator):
        """
        Create a new mathematical expression from the given operator and the last two expressions in self.expression.

        :param operator: Target operator for the new expression.
        :type operator: str
        """
        expr = f"({self.expression.pop()}{operator}{self.expression.pop()})"
        self.expression.append(expr)

    def parse_logical(self, operator):
        """
        Create a new logical expression from the given operator and the last two expressions in self.expression.

        :param operator: Target operator for the new expression.
        :type operator: str
        """
        expr = f"{self.expression.pop()} {operator} {self.expression.pop()}"
        self.expression.append(expr)

    def parse_comparison(self, comparator):
        """
        Create a new comparison expression from the given operator and the last two expressions in self.expression.

        :param operator: Target operator for the new expression.
        :type operator: str
        """
        expr = f"{self.expression.pop()} {comparator} {self.expression.pop()}"
        self.expression.append(expr)

    def visit_Name(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Name nodes.

        :param node: A node of type 'Name' from an AST.
        :type node: ast.Name
        """
        self.expression.append(node.id)
        self.generic_visit(node)

    def visit_Constant(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Constant nodes.

        :param node: A node of type 'Constant' from an AST.
        :type node: ast.Constant
        """
        self.expression.append(str(node.value))
        self.generic_visit(node)

    ###########################################################################
    ### The below methods are deprecated as of Python 3.8.                  ###
    ### They are included for compatibility with 3.7 and earlier.           ###
    ### If <=3.7 becomes unsupported, please remove these visitor methods.  ###
    def visit_Num(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Num nodes.

        :param node: A node of type 'Num' from an AST.
        :type node: ast.Num
        """
        self.expression.append(str(node.n))
        self.generic_visit(node)

    def visit_Str(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Str nodes.

        :param node: A node of type 'Str' from an AST.
        :type node: ast.Str
        """
        self.expression.append(str(node.s))
        self.generic_visit(node)

    def visit_Bytes(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Bytes nodes.

        :param node: A node of type 'Bytes' from an AST.
        :type node: ast.Bytes
        """
        self.expression.append(str(node.s))
        self.generic_visit(node)

    def visit_NameConstant(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.NameConstant nodes.

        :param node: A node of type 'NameConstant' from an AST.
        :type node: ast.NameConstant
        """
        self.expression.append(str(node.value))
        self.generic_visit(node)

    def visit_Ellipsis(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Ellipsis nodes.

        :param node: A node of type 'Ellipsis' from an AST.
        :type node: ast.Ellipsis
        """
        self.expression.append(str(node))
        self.generic_visit(node)
    ### End of deprecated functions (deprecated as of Python 3.8)           ###
    ###########################################################################

    def visit_BinOp(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.BinOp nodes. Right node is visited first.
        By visiting the left node last, the most recently appended token is always the left-hand token.
        This allows us to always append when adding to the expression, and always pop when processing it.

        :param node: A node of type 'BinOp' from an AST.
        :type node: ast.BinOp
        """
        self.visit(node.right)
        self.visit(node.left)
        self.visit(node.op)

    def visit_Call(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Call nodes.

        :param node: A node of type 'Call' from an AST.
        :type node: ast.Call
        """
        arg_list = []
        for arg in node.args:
            self.visit(arg)
            arg_list.append(self.expression.pop())
        arg_list = ",".join(arg_list)
        expr = f"{node.func.id}({arg_list})"
        self.expression.append(expr)

    def visit_BoolOp(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.BoolOp nodes.
        Base converter class assumes that "And" and "Or" operations are defined by inheriting class.
        Implement visit_And() and visit_Or() to define behavior.

        :param node: A node of type 'BoolOp' from an AST.
        :type node: ast.BoolOp
        """
        for operand in reversed(node.values):
            self.visit(operand)
        # Process n-1 operations; for n operands, there are n-1 operations.
        # Example: x && y || z -> 3 operands, 2 operations
        for _ in range(len(node.values) - 1):
            self.visit(node.op)

    def visit_Add(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Add nodes.

        :param node: A node of type 'Add' from an AST.
        :type node: ast.Add
        """
        self.generic_visit(node)
        self.parse_operator("+")

    def visit_Sub(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Sub nodes.

        :param node: A node of type 'Sub' from an AST.
        :type node: ast.Sub
        """
        self.generic_visit(node)
        self.parse_operator("-")

    def visit_Mult(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Mult nodes.

        :param node: A node of type 'Mult' from an AST.
        :type node: ast.Mult
        """
        self.generic_visit(node)
        self.parse_operator("*")

    def visit_Div(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Div nodes.

        :param node: A node of type 'Div' from an AST.
        :type node: ast.Div
        """
        self.generic_visit(node)
        self.parse_operator("/")

    def visit_Pow(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Pow nodes.

        :param node: A node of type 'Pow' from an AST.
        :type node: ast.Pow
        """
        self.generic_visit(node)
        self.parse_operator("**")

    def visit_Compare(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Compare nodes.

        :param node: A node of type 'Compare' from an AST.
        :type node: ast.Compare
        """
        for comparator in node.comparators:
            self.visit(comparator)
        self.visit(node.left)
        for operator in node.ops:
            self.visit(operator)

    def visit_Eq(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Eq nodes.

        :param node: A node of type 'Eq' from an AST.
        :type node: ast.Eq
        """
        self.generic_visit(node)
        self.parse_comparison("==")

    def visit_NotEq(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.NotEq nodes.

        :param node: A node of type 'NotEq' from an AST.
        :type node: ast.NotEq
        """
        self.generic_visit(node)
        self.parse_comparison("!=")

    def visit_Lt(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Lt nodes.

        :param node: A node of type 'Lt' from an AST.
        :type node: ast.Lt
        """
        self.generic_visit(node)
        self.parse_comparison("<")

    def visit_LtE(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.LtE nodes.

        :param node: A node of type 'LtE' from an AST.
        :type node: ast.LtE
        """
        self.generic_visit(node)
        self.parse_comparison("<=")

    def visit_Gt(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Gt nodes.

        :param node: A node of type 'Gt' from an AST.
        :type node: ast.Gt
        """
        self.generic_visit(node)
        self.parse_comparison(">")

    def visit_GtE(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.GtE nodes.

        :param node: A node of type 'GtE' from an AST.
        :type node: ast.GtE
        """
        self.generic_visit(node)
        self.parse_comparison(">=")

    def visit_UnaryOp(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.UnaryOp nodes.

        :param node: A node of type 'UnaryOp' from an AST.
        :type node: ast.UnaryOp
        """
        self.visit(node.operand)
        self.visit(node.op)

    def visit_USub(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.USub nodes.

        :param node: A node of type 'USub' from an AST.
        :type node: ast.USub
        """
        self.expression.append(f"-{self.expression.pop()}")

    def get_str(self):
        """
        Get the string representation of the expression.

        :returns: The string expression of self.tree.
        :rtype: str
        """
        return self._get_str(self.tree)

class PythonConverter(ExpressionConverter):
    """
    Converts an ast.AST to a Python expression string.
    """
    def visit_And(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.And nodes.

        :param node: A node of type 'And' from an AST.
        :type node: ast.And
        """
        self.parse_logical("and")

    def visit_Or(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Or nodes.

        :param node: A node of type 'Or' from an AST.
        :type node: ast.Or
        """
        self.parse_logical("or")

class CppConverter(ExpressionConverter):
    """
    Converts an ast.AST to a C++ expression string.
    """
    class CppExpressionTransformer(ast.NodeTransformer):
        """
        A subclass of ast.NodeTransformer used to convert ast.AST expresions to C++ expressions.
        """
        def visit_BinOp(self, node: "ast.BinOp"): # pylint: disable=invalid-name
            """
            Visitor function for ast.BinOp nodes.

            :param node: A node of type 'BinOp' from an AST.
            :type node: ast.BinOp

            :returns: The original node.
            :rtype: ast.BinOp
            """
            self.generic_visit(node)
            if isinstance(node.op, ast.Pow):
                node = ast.copy_location(ast.Call(
                    func=ast.Name(id='pow', ctx=ast.Load()),
                    args=[node.left, node.right],
                    keywords=[]
                ), node)
            return node

    def visit_And(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.And nodes.

        :param node: A node of type 'And' from an AST.
        :type node: ast.And
        """
        self.generic_visit(node)
        self.parse_logical("&&")

    def visit_Or(self, node): # pylint: disable=invalid-name
        """
        Visitor function for ast.Or nodes.

        :param node: A node of type 'Or' from an AST.
        :type node: ast.Or
        """
        self.generic_visit(node)
        self.parse_logical("||")

    def get_str(self):
        """
        Get the string representation of the expression.

        :returns: The string expression of super().tree.
        :rtype: str
        """
        expr = CppConverter.CppExpressionTransformer().visit(self.tree)
        return super()._get_str(expr)
