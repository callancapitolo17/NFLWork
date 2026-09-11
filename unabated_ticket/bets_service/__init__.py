"""Local bets service for the Unabated Ticket extension (#114).

Serves the user's own bet history as normalised records on
http://127.0.0.1:8094/bets.json so the panel can flag lines already bet.
One `Source` per venue (sources/); Kalshi ships here, #115/#116/#117 add
BetOnline / Novig / ProphetX behind the same protocol.
"""
