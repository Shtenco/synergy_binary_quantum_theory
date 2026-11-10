#!/usr/bin/env python3
"""Independent audit of spectral_mass_theory/mass_relation.json.

Usage:
  python audit_mass_relations.py path/to/mass_relation.json

Only a small arithmetic AST whitelist is evaluated. No eval/exec is used.
"""
from __future__ import annotations
import ast
import json
import math
import sys
from collections import Counter


def safe_eval(expr: str, K: float = 6.0) -> float:
    tree = ast.parse(expr, mode="eval")

    def visit(node):
        if isinstance(node, ast.Expression):
            return visit(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name):
            if node.id == "K":
                return K
            if node.id == "pi":
                return math.pi
            raise ValueError(f"unknown name: {node.id}")
        if isinstance(node, ast.Attribute):
            if isinstance(node.value, ast.Name) and node.value.id == "math" and node.attr == "pi":
                return math.pi
            raise ValueError("unsupported attribute")
        if isinstance(node, ast.Call):
            if (isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "math"
                and node.func.attr == "sqrt"
                and len(node.args) == 1 and not node.keywords):
                return math.sqrt(visit(node.args[0]))
            if isinstance(node.func, ast.Name) and node.func.id == "sqrt" and len(node.args) == 1:
                return math.sqrt(visit(node.args[0]))
            raise ValueError("unsupported function call")
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
            return +visit(node.operand)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add): return a + b
            if isinstance(node.op, ast.Sub): return a - b
            if isinstance(node.op, ast.Mult): return a * b
            if isinstance(node.op, ast.Div): return a / b
            if isinstance(node.op, ast.Pow): return a ** b
        raise ValueError(f"unsupported AST node: {type(node).__name__}")

    return visit(tree)


def main(path: str) -> None:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    particles = data.get("particles", [])
    relations = data.get("relations", [])
    by_name = {p["name"]: p for p in particles}
    by_pair = {r["pair"]: r for r in relations}

    parse_errors = []
    formula_mismatches = []
    teta_mismatches = []
    experimental_ratio_mismatches = []
    reciprocal_pairs = []
    missing_reverse = []
    formula_counts = Counter()

    for row in relations:
        pair = row["pair"]
        formula_counts[row["formula"]] += 1
        try:
            computed = safe_eval(row["formula"], data.get("constants", {}).get("K", 6))
            relerr = abs(computed / row["model_value"] - 1.0)
            if relerr > 1e-10:
                formula_mismatches.append((pair, row["formula"], computed, row["model_value"], relerr))
        except Exception as exc:
            parse_errors.append((pair, row["formula"], repr(exc)))

        # The article prose says experimental/model. Check both conventions.
        model_over_exp = row["model_value"] / row["experimental_value"]
        exp_over_model = row["experimental_value"] / row["model_value"]
        if abs(row["teta"] - exp_over_model) > 1e-10:
            teta_mismatches.append((pair, row["teta"], exp_over_model, model_over_exp))

        names = pair.split(" / ")
        if len(names) == 2 and names[0] in by_name and names[1] in by_name:
            ratio = by_name[names[0]]["mass"] / by_name[names[1]]["mass"]
            if abs(ratio / row["experimental_value"] - 1.0) > 1e-10:
                experimental_ratio_mismatches.append((pair, ratio, row["experimental_value"]))

        if len(names) == 2:
            reverse = by_pair.get(names[1] + " / " + names[0])
            if reverse and names[0] < names[1]:
                reciprocal_pairs.append((pair,
                    row["model_value"] * reverse["model_value"],
                    row["experimental_value"] * reverse["experimental_value"]))
            elif not reverse:
                missing_reverse.append(pair)

    print(f"particles={len(particles)} unique_particle_names={len(by_name)}")
    print(f"relations={len(relations)} expected_directed_pairs={len(particles)*(len(particles)-1)}")
    print(f"unique_formula_strings={len(formula_counts)}")
    print(f"parse_errors={len(parse_errors)}")
    print(f"literal_formula_value_mismatches={len(formula_mismatches)}")
    print(f"teta_mismatches_vs_article_definition={len(teta_mismatches)}")
    print(f"experimental_ratio_mismatches_vs_particle_mass_records={len(experimental_ratio_mismatches)}")
    print(f"reciprocal_unordered_pairs={len(reciprocal_pairs)} missing_reverse={len(missing_reverse)}")
    if reciprocal_pairs:
        print("max_model_reciprocal_product_error=",
              max(abs(x[1] - 1.0) for x in reciprocal_pairs))
        print("max_experimental_reciprocal_product_error=",
              max(abs(x[2] - 1.0) for x in reciprocal_pairs))
    for title, rows in [
        ("PARSE ERRORS", parse_errors),
        ("FORMULA MISMATCHES", sorted(formula_mismatches, key=lambda x: -x[4])[:12]),
        ("TETA MISMATCHES", teta_mismatches[:5]),
        ("EXPERIMENTAL RATIO MISMATCHES", experimental_ratio_mismatches[:5]),
    ]:
        print("\n" + title)
        for row in rows:
            print(row)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python audit_mass_relations.py path/to/mass_relation.json")
    main(sys.argv[1])
