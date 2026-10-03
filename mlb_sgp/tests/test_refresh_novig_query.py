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


def _bundle(map_sources=(QUERY_SRC, FRAG_B_SRC, FRAG_A_SRC), ast=AST_LITERAL,
            raw_entries=()):
    names = {QUERY_SRC: "EventMarkets_QueryDocument",
             FRAG_A_SRC: "A_FragFragmentDoc", FRAG_B_SRC: "B_FragFragmentDoc"}
    entries = [f"{json.dumps(src)}:_.{names[src]}" for src in map_sources]
    entries.extend(raw_entries)
    return f"var docs={{{','.join(entries)}}};{ast};"


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


def test_an_undecodable_unrelated_document_does_not_block_the_refresh():
    """The minifier writes Latin-1 characters as JS-only \\xNN escapes,
    which JSON rejects; one in a document we don't need is skipped."""
    other = r'"\n  query Other_Query {\n    caf\xe9\n  }\n":_.Other_QueryDocument'
    bundle = _bundle(raw_entries=(other,))
    assert refresh.build_query(bundle) == refresh.build_query(_bundle())


def test_declared_variables_reads_the_operation_header():
    query = ("query EventMarkets_Query($eventId: uuid, $marketVisibleWhere: "
             "market_bool_exp) @cached(ttl: 5) {\n  event(where: {id: "
             "{_eq: $eventId}}) { id }\n}")
    assert refresh.declared_variables(query) == {"eventId", "marketVisibleWhere"}


def _leg(leg_id, available=0.5):
    return {"id": leg_id, "available": available}


def _structure(**overrides):
    structure = {"home_spread": {-1.5: _leg("hs")}, "away_spread": {1.5: _leg("as")},
                 "over": {8.5: _leg("ov")}, "under": {8.5: _leg("un")},
                 "home_ml": None, "away_ml": None}
    structure.update(overrides)
    return structure


def test_two_sided_rungs_accepts_a_priceable_tree():
    assert refresh.has_two_sided_rungs(_structure())


def test_two_sided_rungs_rejects_a_tree_missing_prices_or_a_side():
    no_prices = _structure(over={8.5: _leg("ov", available=None)})
    one_sided_spread = _structure(away_spread={})
    assert not refresh.has_two_sided_rungs(no_prices)
    assert not refresh.has_two_sided_rungs(one_sided_spread)
