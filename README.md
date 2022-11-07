# WaterSortCpp

This is a simple yet very fast solver for the Water Sort puzzle. It uses the a* algorithm to in order to find a solution. 
The heuristic function it uses in order to find the distance from the solved state is not perfect, therefore in some rare cases the solution found is not optimal - the one with the fewer steps.

## Example usage

```
main()
{
	// Initialize a new game with 12 colors, 2 empty bottles
	// Height of each bottle is 4
	game game(12, 2, 4);

	// Generate an initial random state
	state st = game.random_state(123);
	
  	// Initialize the parent node
	node initial(st, nullptr);

	// Find the solution
	node sol = game.solve(initial);

	// Print all the steps to solution
	cout << game.get_history(sol) << endl;
}

```
Expected output:
```
|4421|721B|4632|140B|5608|9965|9623|B0A3|A590|7387|B788|A1A5|    |    |  0 38
|4421|721B|463 |140B|5608|9965|9623|B0A3|A590|7387|B788|A1A5|2   |    |  1 37
|4421|721B|4633|140B|5608|9965|962 |B0A3|A590|7387|B788|A1A5|2   |    |  2 36
|4421|721B|4633|140B|5608|9965|96  |B0A3|A590|7387|B788|A1A5|22  |    |  3 35
|4421|721B|4633|140B|560 |9965|96  |B0A3|A590|7387|B788|A1A5|22  |8   |  4 34
|4421|721B|4633|140B|560 |9965|96  |B0A3|A590|7387|B7  |A1A5|22  |888 |  5 33
|4421|721B|4633|140B|560 |9965|96  |B0A3|A590|738 |B77 |A1A5|22  |888 |  6 32
|4421|721B|4633|140B|560 |9965|96  |B0A3|A590|73  |B77 |A1A5|22  |8888|  7 31
|4421|721B|46  |140B|560 |9965|96  |B0A3|A590|7333|B77 |A1A5|22  |8888|  8 30
|4421|721B|46  |140B|5600|9965|96  |B0A3|A59 |7333|B77 |A1A5|22  |8888|  9 29
|4421|721B|466 |140B|5600|9965|9   |B0A3|A59 |7333|B77 |A1A5|22  |8888| 10 28
|4421|721B|466 |140B|5600|9965|    |B0A3|A599|7333|B77 |A1A5|22  |8888| 11 27
|4421|721B|466 |140B|5600|9965|3   |B0A |A599|7333|B77 |A1A5|22  |8888| 12 26
|4421|721B|466 |140B|5600|9965|3333|B0A |A599|7   |B77 |A1A5|22  |8888| 13 25
|4421|721B|466 |140B|5600|9965|3333|B0A |A599|    |B777|A1A5|22  |8888| 14 24
|4421|721B|466 |140B|56  |9965|3333|B0A |A599|00  |B777|A1A5|22  |8888| 15 23
|4421|721B|4666|140B|5   |9965|3333|B0A |A599|00  |B777|A1A5|22  |8888| 16 22
|4421|721B|4666|140B|55  |996 |3333|B0A |A599|00  |B777|A1A5|22  |8888| 17 21
|4421|721B|4666|140B|555 |996 |3333|B0A |A599|00  |B777|A1A |22  |8888| 18 20
|4421|721B|4666|140B|555 |996 |3333|B0  |A599|00  |B777|A1AA|22  |8888| 19 19
|4421|721B|4666|140B|555 |996 |3333|B   |A599|000 |B777|A1AA|22  |8888| 20 18
|4421|721B|4666|140 |555 |996 |3333|BB  |A599|000 |B777|A1AA|22  |8888| 21 17
|4421|721 |4666|140 |555 |996 |3333|BBB |A599|000 |B777|A1AA|22  |8888| 22 16
|442 |7211|4666|140 |555 |996 |3333|BBB |A599|000 |B777|A1AA|22  |8888| 23 15
|44  |7211|4666|140 |555 |996 |3333|BBB |A599|000 |B777|A1AA|222 |8888| 24 14
|44  |7211|4666|14  |555 |996 |3333|BBB |A599|0000|B777|A1AA|222 |8888| 25 13
|    |7211|4666|1444|555 |996 |3333|BBB |A599|0000|B777|A1AA|222 |8888| 26 12
|666 |7211|4   |1444|555 |996 |3333|BBB |A599|0000|B777|A1AA|222 |8888| 27 11
|666 |7211|4444|1   |555 |996 |3333|BBB |A599|0000|B777|A1AA|222 |8888| 28 10
|6666|7211|4444|1   |555 |99  |3333|BBB |A599|0000|B777|A1AA|222 |8888| 29  9
|6666|7211|4444|1   |555 |9999|3333|BBB |A5  |0000|B777|A1AA|222 |8888| 30  8
|6666|7211|4444|1   |5555|9999|3333|BBB |A   |0000|B777|A1AA|222 |8888| 31  7
|6666|7211|4444|1   |5555|9999|3333|BBB |AAA |0000|B777|A1  |222 |8888| 32  6
|6666|72  |4444|1   |5555|9999|3333|BBB |AAA |0000|B777|A111|222 |8888| 33  5
|6666|7   |4444|1   |5555|9999|3333|BBB |AAA |0000|B777|A111|2222|8888| 34  4
|6666|7777|4444|1   |5555|9999|3333|BBB |AAA |0000|B   |A111|2222|8888| 35  3
|6666|7777|4444|1   |5555|9999|3333|BBBB|AAA |0000|    |A111|2222|8888| 36  2
|6666|7777|4444|1111|5555|9999|3333|BBBB|AAA |0000|    |A   |2222|8888| 37  1
|6666|7777|4444|1111|5555|9999|3333|BBBB|AAAA|0000|    |    |2222|8888| 38  0
```
An alphabet is used to represent colors. When a state is created, it is sorted using this alphabet. This helps eliminate equivalent states, that is, states that differ only by the order of the bottles.


