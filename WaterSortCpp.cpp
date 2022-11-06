#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include <sstream>
#include <iostream>
#include <chrono>

using namespace std;

typedef vector<char> bottle;
typedef vector<bottle> state;

/// <summary>
/// Represents a state together with its generation, its heuristing distance 
/// to the final solution, plus a reference to its parent node.
/// </summary>
class node
{

public:
	state st;	// The state of the node (bottles and their content)
	int g;		// Generation
	int h;		// Heuristic estimation of steps to solution
	const node* parent;

	/// <summary>
	/// Constructor of the new node
	/// </summary>
	node(state state, const node* p) : st(state), parent(p)
	{
		g = (p == nullptr) ? 0 : p->g + 1;
		compute_h();
	}

	int f() const
	{
		return g + h;
	}

	/// <summary>
	/// Computation of heuristic distance to solution
	/// </summary>
	void compute_h()
	{
		h = 0;
		for (int i = 0; i < st.size(); i++)
		{
			if (st[i].empty()) continue;

			for (int j = 0; j < st[i].size() - 1; j++)
			{
				if (st[i][j] != st[i][j + 1])
					h++;
			}
			for (int k = i + 1; k < st.size(); k++)
			{
				if (st[k].size() > 0 && st[i][0] == st[k][0])
					h++;
			}
		}
	}
};


/// <summary>
/// Node comparer. Comparison is by their f value, then by distance and then by their states
/// </summary>
/// <param name="a">First node to compare</param>
/// <param name="b">Second node to compare</param>
/// <returns>Returns true if node a is "less than" node b</returns>
struct ncomp
{
	bool operator() (const node& a, const node& b) const
	{
		int fa = a.g + a.h;
		int fb = b.g + b.h;

		return (fa == fb) ? ((a.g == b.g) ? a.st < b.st : a.h < b.h) : (fa < fb);

		//// this gives slightly less steps to solve but it is considerably slower in execution
		//return (fa == fb) ? ((a.g == b.g) ? a.st < b.st : a.g < b.g) : (fa < fb);
	}
};


/// <summary>
/// Represents a game setting, e.g. number of colors, height of each bottle,
/// number of empty bottles
/// </summary>
class game
{
	static string alphabet; 

	int _t;	// Number of total bottles
	int _n; // Number of filled bottles
	int _k; // Number of empty bottles
	int _h; // Height of each bottle

public:
	set<node, ncomp> open;
	set<node, ncomp> closed;

	/// <summary>
	/// Generates a new game 
	/// </summary>
	/// <param name="n">Number of colors</param>
	/// <param name="k">Number of empty bottles</param>
	/// <param name="h">Height of each bottle</param>
	game(int n, int k, int h) : _t(n + k), _n(n), _h(h), _k(k)
	{
	}

	/// <summary>
	/// Generates randomly a new state using a seed value.
	/// </summary>
	state random_state(unsigned int seed)
	{
		vector<char> chars;
		for (int i = 0; i < _n; i++)
		{
			for (int j = 0; j < _h; j++)
			{
				chars.push_back(alphabet[i]);
			}
		}
		auto rng = std::default_random_engine{ seed };
		std::shuffle(chars.begin(), chars.end(), rng);
		int idx = 0;
		state state;
		for (int i = 0; i < _n; i++)
		{
			bottle b;
			for (int j = 0; j < _h; j++)
			{
				b.push_back(chars[idx++]);
			}
			state.push_back(b);
		}
		for (int i = 0; i < _k; i++)
		{
			state.push_back(bottle());
		}
		sort(state.begin(), state.end());
		return state;
	}

	/// <summary>
	/// Generates the child nodes of a node
	/// </summary>
	vector<node> get_child_nodes(const node& n)
	{
		vector<node> result;
		const state& st = n.st;
		for (int from = 0; from < _t; from++)
		{
			if (st[from].size() == 0) continue;

			for (int to = 0; to < _t; to++)
			{
				if (to == from || st[to].size() == _h)
					continue;

				state nst(st);
				bool changed = false;
				while (nst[from].size() > 0
					&& nst[to].size() < _h
					&& (nst[to].size() == 0 ||
						(nst[to].back() == nst[from].back()))
					)
				{
					changed = true;
					nst[to].push_back(nst[from].back());
					nst[from].pop_back();
				}
				if (changed)
				{
					sort(nst.begin(), nst.end());
					result.push_back(node(nst, &n));
				}
			}
		}
		return result;
	}

	/// <summary>
	/// Checks whether we can skip entering a node in the 
	/// open or in the closed set
	/// </summary>
	/// <param name="start">The set containing other nodes</param>
	/// <returns>The node to check</returns>
	bool check_skip(const set<node, ncomp>& s, const node& n)
	{
		auto it = s.begin();
		while (it != s.end() && (*it).f() < n.f())
		{
			if ((*it).st == n.st) return true;
			it++;
		}
		return false;
	}

	/// <summary>
	/// Returns a string representation of a node
	/// </summary>
	/// <param name="start">The node to convert to a string</param>
	/// <returns>A string for the node state, generation and heuristic distance</returns>
	string node_to_string(const node& n)
	{
		ostringstream ostr;
		ostr << '|';
		for (int i = 0; i < _t; i++)
		{
			for (int j = 0; j < _h; j++)
			{
				if (j < n.st[i].size())
					ostr << n.st[i][j];
				else
					ostr << " ";
			}
			ostr << '|';
		}
		ostr.width(3);
		ostr << n.g;
		ostr.width(3);
		ostr << n.h;
		return ostr.str();
	}

	/// <summary>
	/// Builds recursively a string representation of a node 
	/// together with all its parent nodes
	/// </summary>
	string get_history(const node& n)
	{
		if (n.parent == nullptr)
		{
			return node_to_string(n);
		}
		else
		{
			return get_history(*n.parent) + "\n" + node_to_string(n);
		}
	}

	/// <summary>
	/// Solves a Water Sort problem 
	/// </summary>
	node solve(node start)
	{
		open.insert(start);
		while (open.size() > 0)
		{
			auto it = open.begin();
			node n = *it;
			open.erase(it);

			// we have reached a solution!
			if (n.h == 0)
				return n;

			auto p = closed.insert(n);
			if (p.second)
			{
				auto& parent = *p.first;
				vector<node> children = get_child_nodes(parent);

				for (node& item : children)
				{
					if (check_skip(open, item))
						continue;
					if (check_skip(closed, item))
						continue;
					open.insert(item);
				}
			}
		}
		// return a node with an empty state to indicate
		// that the problem is unsolved
		return node(state(), nullptr);
	}
};

// Extend this if you plan to use more than 16 colors
string game::alphabet = "0123456789ABCDEF";

int main()
{
	// Initialize a new game with 12 colors, 2 empty bottles
	// Height of each bottle is 4
	game game(12, 2, 4);

	// Generate an initial random state
	state st = game.random_state(123);
	// Initialize parent node
	node initial(st, nullptr);

	// Solve the state and get the solution
	node sol = game.solve(initial);

	// Print the steps 
	cout << game.get_history(sol) << endl;

	// Uncomment the following section in order to run
	// the solver several times

	//int N = 100;
	//int c = 0;
	//double sum = 0.0;
	//auto start = chrono::high_resolution_clock::now();
	//for (int i = 0; i < N; i++)
	//{
	//	game g(12, 2, 4);
	//	state st = g.random_state(i);
	//	node n(st, nullptr);
	//	node solution = g.solve(n);
	//	if (solution.st.size() > 0)
	//	{
	//		cout.width(3);
	//		cout << i;
	//		cout.width(4);
	//		cout << solution.g  << endl;
	//		sum += solution.g;
	//		c++;
	//	}
	//	else
	//	{
	//		cout.width(3);
	//		cout << i << "   Unsolved!" << endl;
	//	}
	//}
	//auto end = chrono::high_resolution_clock::now();
	//chrono::duration<double> d = end - start;
	//std::cout << "Average steps = " << sum / c << "  Average seconds = " << d.count() / c << endl;
}
