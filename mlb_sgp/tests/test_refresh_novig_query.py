"""Unit tests for refresh_novig_query's bundle extraction — no network.

The synthetic bundle mimics the two pieces graphql-codegen's client preset
puts in the novig.com bundle: a documents map keyed by each definition's
JS-escaped source text, and the operation's pre-parsed AST literal whose
`definitions` fix the order Novig's allowlist requires.
"""
import json

import pytest

from mlb_sgp import refresh_novig_query as refresh

QUERY_SRC = ("\n  query EventMarkets_Query($eventId: uuid) {\n"
             "    event(where: { id: { _eq: $eventId } }) {\n"
             "      ...B_Frag\n    }\n  }\n")
FRAG_A_SRC = "\n  fragment A_Frag on event {\n    id\n  }\n"
FRAG_B_SRC = "\n  fragment B_Frag on event {\n    ...A_Frag\n    type\n  }\n"

# AST order (operation, A, B) differs from the order the map lists them in
# (operation, B, A), so a test passing proves the AST decides. The StringValue
# with braces in it checks the bracket matcher skips string contents.
AST_LITERAL = (
    'ht={kind:"Document",definitions:['
    '{kind:"OperationDefinition",operation:"query",'
    'name:{kind:"Name",value:"EventMarkets_Query"},'
    'directives:[{kind:"Directive",arguments:[{kind:"Argument",'
    'value:{kind:"StringValue",value:"a}]b",block:!1}}]}]},'
    '{kind:"FragmentDefinition",name:{kind:"Name",value:"A_Frag"},'
    'selectionSet:{kind:"SelectionSet",selections:[]}},'
    '{kind:"FragmentDefinition",name:{kind:"Name",value:"B_Frag"},'
    'selectionSet:{kind:"SelectionSet",selections:[]}}]},'
    'Tt={kind:"Document",definitions:[{kind:"FragmentDefinition",'
    'name:{kind:"Name",value:"Unrelated_Frag"}}]}'
)


def _bundle(map_sources=(QUERY_SRC, FRAG_B_SRC, FRAG_A_SRC), ast=AST_LITERAL):
    names = {QUERY_SRC: "EventMarkets_QueryDocument",
             FRAG_A_SRC: "A_FragFragmentDoc", FRAG_B_SRC: "B_FragFragmentDoc"}
    entries = ",".join(f"{json.dumps(src)}:_.{names[src]}" for src in map_sources)
    return f"var docs={{{entries}}};{ast};"


def test_definition_order_comes_from_the_ast_literal():
    assert refresh.definition_order(_bundle()) == [
        "EventMarkets_Query", "A_Frag", "B_Frag"]


def test_build_query_joins_sources_in_ast_order():
    expected = "\n\n".join(s.strip() for s in (QUERY_SRC, FRAG_A_SRC, FRAG_B_SRC))
    assert refresh.build_query(_bundle()) == expected


def test_build_query_fails_when_a_fragment_source_is_missing():
    with pytest.raises(RuntimeError, match="missing sources=\\['A_Frag'\\]"):
        refresh.build_query(_bundle(map_sources=(QUERY_SRC, FRAG_B_SRC)))


def test_definition_order_fails_without_the_ast_literal():
    with pytest.raises(RuntimeError, match="found 0"):
        refresh.definition_order(_bundle(ast="ht={}"))
