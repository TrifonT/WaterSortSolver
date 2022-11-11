#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include <sstream>
#include <iostream>
#include <chrono>
#include <map>

using namespace std;

class bottle : public vector<char>
{
public:
	int order = 0; // position of the bottle in a state
};

typedef vector<bottle> state;

/// <summary>
/// Represents a state together with its generation, its heuristing distance 
/// to the final solution, plus a reference to its parent node.
/// </summary>
class node
{
public:
	state st;	// The state of the node (bottles and their content)
	mutable int g;		// Generation
	int h;		// Heuristic estimation of number of steps to reach solution
	const node* parent;

	/// <summary>
	/// Constructor
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
	/// Computation of a heuristic distance to the solution
	/// (number of steps needed) 
	/// </summary>
	void compute_h()
	{
		static map<char, int> m;
		m.clear();
		h = 0;
		for (size_t i = 0; i < st.size(); i++)
		{
			if (st[i].empty()) continue;

			for (size_t j = 0; j < st[i].size() - 1; j++)
			{
				if (st[i][j] != st[i][j + 1])
					h++;
			}
			m[st[i][0]]++;
		}
		for (auto& item : m)
		{
			h += item.second - 1;
		}
	}
};

/// <summary>
/// Node comparer for the open set. Compares nodes by their f value,
/// then by distance and then by their states
/// </summary>
/// <param name="a">First node to compare</param>
/// <param name="b">Second node to compare</param>
/// <returns>Returns true if node a is "less than" node b</returns>
struct ocomp
{
	bool operator() (const node& a, const node& b) const
	{
		int fa = a.g + a.h;
		int fb = b.g + b.h;

		return (fa == fb) ? ((a.g == b.g) ? a.st < b.st : a.h < b.h) : (fa < fb);
	}
};

/// <summary>
/// Node comparer for the closed set. Compares nodes by their states
/// only
/// </summary>
/// <param name="a">First node to compare</param>
/// <param name="b">Second node to compare</param>
/// <returns>Returns true if node a is "less than" node b</returns>
struct ccomp
{
	bool operator() (const node& a, const node& b) const
	{
		return a.st < b.st;
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
	set<node, ocomp> open;
	set<node, ccomp> closed;

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
			b.order = i;
			for (int j = 0; j < _h; j++)
			{
				b.push_back(chars[idx++]);
			}
			state.push_back(b);
		}
		for (int i = 0; i < _k; i++)
		{
			bottle b;
			b.order = _n + i;
			state.push_back(b);
		}
		sort(state.begin(), state.end());
		return state;
	}

	/// <summary>
	/// Returns a string representation of a node
	/// </summary>
	/// <param name="start">The node to convert to a string</param>
	/// <returns>A string for the node state, generation and heuristic distance</returns>
	string node_to_string(const node& n)
	{
		state stc(n.st);
		// We sort bottles according to their initial order
		sort(stc.begin(), stc.end(), [&](bottle a, bottle b) {return a.order < b.order; });

		ostringstream ostr;
		ostr << '|';
		for (int i = 0; i < _t; i++)
		{
			for (int j = 0; j < _h; j++)
			{
				if (j < stc[i].size())
					ostr << stc[i][j];
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
	/// together with all its anchestor nodes
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
	/// Generates the child nodes of a node
	/// </summary>
	vector<node> get_child_nodes(const node& n)
	{
		vector<node> result;		
		const state& st = n.st;
		for (size_t from = 0; from < _t; from++)
		{
			if (st[from].empty()) continue;

			for (size_t to = 0; to < _t; to++)
			{
				if (to == from || st[to].size() == _h)
					continue;

				if (st[to].empty() || st[from].back() == st[to].back())
				{					
					state nst(st);
					do
					{
						nst[to].push_back(nst[from].back());
						nst[from].pop_back();
					} while (
						(!nst[from].empty()) &&
						(nst[to].size() < _h) &&
						(nst[to].back() == nst[from].back())
						);
					sort(nst.begin(), nst.end());
					result.emplace_back(nst, &n);					
				}
			}
		}
		return result;
	}

	/// <summary>
	/// Checks whether we can skip entering a node in the 
	/// open set
	/// </summary>
	/// <param name="start">The set containing other nodes</param>
	/// <returns>The node to check</returns>
	bool check_skip(const set<node, ocomp>& s, const node& n)
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
	/// Checks whether we can skip entering a node in the 
	/// closed set
	/// </summary>
	/// <param name="start">The set containing other nodes</param>
	/// <returns>The node to check</returns>
	bool check_skip(const set<node, ccomp>& s, const node& n)
	{
		auto it = s.find(n);
		return  (it == s.end()) ? false : n.g >= (*it).g;
	}

	/// <summary>
	/// Solves a Water Sort problem using a* algorithm
	/// </summary>
	node solve(node start)
	{
		// put initial node in open set
		open.insert(start);
		while (open.size() > 0)
		{
			// nodes in the open set are sorted by their f=g+h value
			// here we get the node with the lower f value
			auto it = open.begin();
			node n = *it;
			open.erase(it);

			// we have reached a solution!
			if (n.h == 0)
				return n;

			// insert the node in the closed set
			auto p = closed.emplace(n);

			// if insertion was successful
			if (p.second)
			{
				// get a reference of the node in its new location
				auto& parent = *p.first;
				vector<node> children = get_child_nodes(parent);

				// insert each child node in the open set 
				// if inclusion checks are met
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
		// return a node with an empty state in order to indicate
		return node(state(), nullptr);
	}

};

// Extend this if you plan to use more than 25 colors
string game::alphabet = "0123456789ABCDEFGHIJKLMNO";


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
	//	the solver several timesand collect statistics

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
	//	if (!solution.st.empty())
	//	{
	//		cout.width(3);
	//		cout << i;
	//		cout.width(4);
	//		cout << solution.g << endl;
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
	//std::cout << "Average steps = " << sum / c << "  Average seconds to solve = " << d.count() / N << endl;
}

