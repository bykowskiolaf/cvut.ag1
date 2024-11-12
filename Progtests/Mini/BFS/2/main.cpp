/*
 * # State space exploration using BFS
 *
 * Your task is to implement (template of) the function `solve`
 * which uses BFS to explore a state space. The problem statement
 * typically doesn't explicitly describe the possible states,
 * but it does specify which actions can be performed in each
 * state. Performing an action leads to a new state, and the
 * sequence of performed actions is usually more interesting
 * than the list of visited states. The argument of `solve`
 * function is an object describing the state space. This object
 * contains (see `LostHero` which implements this interface):
 *
 * - Type `State` which represents a possible state (i.e.,
 *   a vertex).
 * - Type `Action` which represents possible action.
 *
 * - Method `initial_state()` which return the initial state.
 *
 * - Method `is_target(state)` which checks if the given state
 *   is a target state.
 *
 * - Method `next_state(cur_state, action)` which calculates
 *   state obtained by performing action `action` in `cur_state`.
 *
 * - Method `possible_actions(state)` which returns a list of
 *   all actions possible in the state `state`.
 *
 * It returns the sortest possible sequence of actions to reach
 * a target state from the initial state.
 *
 * Then you should solve three problems by modeling their state
 * spaces:
 *
 * - SuperstitiousHero
 * - LimpingKnight
 * - BottleProblem
 *
 * See the template for their description. Do not forget to
 * change macro `TEST_ProblemName` to value `1` after implementing
 * given class.
 *
 * The time limit is 10 seconds for each test.
 *
 */

#ifndef __PROGTEST__
#include <cassert>
#include <cstdarg>
#include <iomanip>
#include <cstdint>
#include <iostream>
#include <memory>
#include <limits>
#include <optional>
#include <array>
#include <algorithm>
#include <deque>
#include <queue>
#include <random>
#include <ranges>
#include <type_traits>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <compare>
#include <ranges>
#include <variant>


struct TestFailed : std::runtime_error {
  using std::runtime_error::runtime_error;
};

std::string fmt(const char *f, ...) {
  va_list args1;
  va_list args2;
  va_start(args1, f);
  va_copy(args2, args1);

  std::string buf(vsnprintf(nullptr, 0, f, args1), '\0');
  va_end(args1);

  vsnprintf(buf.data(), buf.size() + 1, f, args2);
  va_end(args2);

  return buf;
}

#define CHECK(succ, ...) do { \
    if (!(succ)) throw TestFailed(fmt(__VA_ARGS__)); \
  } while (0)


template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };


constexpr size_t NO_SOLUTION = -size_t(1);


// For BottleProblem
using Bottle = size_t;
struct Fill { Bottle index; };
struct Empty { Bottle index; };
struct Transfer {
  Bottle from;
  Bottle to;
};

using BottleOp = std::variant<Fill, Empty, Transfer>;


// For maze problems
enum struct Tile : uint8_t {
  EMPTY = 0,
  WALL,
  EXIT,
};

using Maze = std::vector<std::vector<Tile>>;

struct Position {
  Position() : Position(-1, -1) {}
  Position(size_t r, size_t c) : row(r), col(c) {}
  friend auto operator <=> (Position, Position) = default;

  size_t row, col;
};

enum struct Direction { UP, DOWN, LEFT, RIGHT };

Position perform_move(Position p, Direction d) {
  auto [ r, c ] = p;
  using enum Direction;
  switch (d) {
    case UP:    return { r - 1, c     };
    case DOWN:  return { r + 1, c     };
    case LEFT:  return { r    , c - 1 };
    case RIGHT: return { r    , c + 1 };
  }
  assert(false);
}

enum struct LimpingKnightAction {
  REST,
  UP, DOWN, LEFT, RIGHT,
  UP_UP_LEFT,
  UP_UP_RIGHT,
  DOWN_DOWN_LEFT,
  DOWN_DOWN_RIGHT,
  LEFT_LEFT_UP,
  LEFT_LEFT_DOWN,
  RIGHT_RIGHT_UP,
  RIGHT_RIGHT_DOWN,
};

bool is_knight_move(LimpingKnightAction a) {
  using enum LimpingKnightAction;
  switch (a) {
    case REST:
    case UP:
    case DOWN:
    case LEFT:
    case RIGHT:
      return false;
    default:
      return true;
  }
}

Position perform_move(Position p, LimpingKnightAction d) {
  auto [ r, c ] = p;
  using enum LimpingKnightAction;
  switch (d) {
    case REST:             return { r    , c     };
    case UP:               return { r - 1, c     };
    case DOWN:             return { r + 1, c     };
    case LEFT:             return { r    , c - 1 };
    case RIGHT:            return { r    , c + 1 };
    case UP_UP_LEFT:       return { r - 2, c - 1 };
    case UP_UP_RIGHT:      return { r - 2, c + 1 };
    case DOWN_DOWN_LEFT:   return { r + 2, c - 1 };
    case DOWN_DOWN_RIGHT:  return { r + 2, c + 1 };
    case LEFT_LEFT_UP:     return { r - 1, c - 2 };
    case LEFT_LEFT_DOWN:   return { r + 1, c - 2 };
    case RIGHT_RIGHT_UP:   return { r - 1, c + 2 };
    case RIGHT_RIGHT_DOWN: return { r + 1, c + 2 };
  }
  assert(false);
}

#endif


// Change to 1 to enable given test
// Change to 1 to enable given test
#define TEST_SuperstitiousHero 1
#define TEST_LimpingKnight 1
#define TEST_BottleProblem 1

// Simple hero in a maze problem
struct LostHero {
  // State is just the position of the hero.
  using State = Position;
  using Action = Direction;

  LostHero(Maze maze, Position hero)
    : _maze(std::move(maze)), _hero(hero) {}

  State initial_state() const { return _hero; }

  // Returns state obtained from v by action a.
  // This method may assume that action a is allowed in state v
  State next_state(State v, Action a) const {
    return perform_move(v, a);
  }

  // Returns list of possible actions in given state
  std::vector<Action> possible_actions(State v) const {
    std::vector<Action> ret;

    using enum Direction;
    for (Direction d : { DOWN, UP, RIGHT, LEFT }) {
      auto [ r, c ] = next_state(v, d);
      if (r >= _maze.size() ||
          c >= _maze[r].size() ||
          _maze[r][c] == Tile::WALL) continue;

      ret.push_back(d);
    }

    return ret;
  }

  bool is_target(State v) const { return _maze[v.row][v.col] == Tile::EXIT; }

  private:
  Maze _maze;
  Position _hero;
};


// Implementation of BFS to solve state space search kind of problems.
// Because the state space is usually much bigger than length of the
// solution, we should not allocate P and D as vectors but we should
// use std::set instead.
//
// We are not interested in vertices from initial state to target one
// but in actions that created this path. The return type is an optional
// to distinguish "no solution exists" (return empty optional) and
// "initial state is a solution" (return optional with empty vector)
// cases.
template<typename SearchProblem>
std::optional<std::vector<typename SearchProblem::Action> > solve(const SearchProblem &G) {
  using State = typename SearchProblem::State;
  using Action = typename SearchProblem::Action;

  std::queue<State> queue;
  std::set<State> visited;
  std::map<State, std::pair<State, Action> > parent;

  State initial = G.initial_state();

  if (G.is_target(initial)) {
    return std::vector<Action>{}; // Empty path
  }

  queue.push(initial);
  visited.insert(initial);

  while (!queue.empty()) {
    State current = queue.front();
    queue.pop();

    for (const Action &action: G.possible_actions(current)) {
      State next = G.next_state(current, action);
      if (visited.count(next)) {
        continue;
      }
      parent[next] = {current, action};
      if (G.is_target(next)) {
        // Reconstruct path
        std::vector<Action> path;
        State s = next;
        while (s != initial) {
          // Get action that led to state s and its parent
          auto [prev_state, act] = parent[s];
          path.push_back(act);
          s = prev_state;
        }
        // Reverse path to get it from initial to target
        std::reverse(path.begin(), path.end());
        return path;
      }
      visited.insert(next);
      queue.push(next);
    }
  }

  return {}; // No solution
}


// Similar problem to LostHero
//
// Moves and maze are the same but the hero is very
// superstitious so his every 5th move must be UP
// (counting moves from 1).
//
// - How will this condition change the definition of State?
// - Is still the shortest way out of the maze a path in the
//   maze in the graph theoretic sense?
#if TEST_SuperstitiousHero
struct SuperstitiousHero {
  struct State {
    Position pos;
    int move_count_mod5{};  // 0 to 4
    friend auto operator<=>(const State&, const State&) = default;
  };
  using Action = Direction;

  SuperstitiousHero(Maze maze, Position hero)
      : _maze(std::move(maze)), _hero(hero) {}

  State initial_state() const { return State{_hero, 0}; }

  std::vector<Action> possible_actions(State v) const {
    std::vector<Action> ret;

    const int next_move_num = v.move_count_mod5 + 1;
    bool must_move_up = (next_move_num % 5 == 0);

    if (must_move_up) {
      // Only consider UP
      Direction d = Direction::UP;
      Position next_pos = perform_move(v.pos, d);
      size_t r = next_pos.row;
      size_t c = next_pos.col;
      if (r < _maze.size() && c < _maze[r].size() && _maze[r][c] != Tile::WALL) {
        ret.push_back(d);
      }
    } else {
      // Consider all possible moves as in LostHero
      using enum Direction;
      for (Direction d : {DOWN, UP, RIGHT, LEFT}) {
        Position next_pos = perform_move(v.pos, d);
        size_t r = next_pos.row;
        size_t c = next_pos.col;
        if (r >= _maze.size() || c >= _maze[r].size() || _maze[r][c] == Tile::WALL)
          continue;
        ret.push_back(d);
      }
    }

    return ret;
  }

  State next_state(State v, const Action& a) const {
    Position next_pos = perform_move(v.pos, a);
    int next_move_count_mod5 = (v.move_count_mod5 + 1) % 5;
    return State{next_pos, next_move_count_mod5};
  }

  bool is_target(State v) const {
    return _maze[v.pos.row][v.pos.col] == Tile::EXIT;
  }

private:
  Maze _maze;
  Position _hero;
};
#endif


// Hero in a maze again
//
// This time the hero is allowed to move as a chess knight (this move
// may skip over walls) but at most once every 8 moves. Other moves
// must be to adjacent tiles as in the previous cases.
// He is also allowed rest (i.e., not move so it is possible
// to make a knight move, rest 7 turns and make a knight move again).
#if TEST_LimpingKnight
struct LimpingKnight {
    struct State {
        Position pos;
        int knight_cooldown; // Values from 0 to 7
        friend auto operator<=>(const State &, const State &) = default;
    };

    using Action = LimpingKnightAction;

    LimpingKnight(Maze maze, Position hero)
        : _maze(std::move(maze)), _hero(hero) {
    }

    // Initial state with position and cooldown set to 0
    State initial_state() const {
        return {_hero, 0};
    }

    // Return all possible actions in the current state
    std::vector<Action> possible_actions(State v) const {
        std::vector<Action> ret;

        // Add adjacent moves and REST
        using enum LimpingKnightAction;
        for (LimpingKnightAction a: {
                 REST,
                 UP,
                 DOWN,
                 LEFT,
                 RIGHT
             }) {
            Position next_pos = perform_move(v.pos, a);
            if (next_pos.row >= _maze.size() ||
                next_pos.col >= _maze[next_pos.row].size() ||
                _maze[next_pos.row][next_pos.col] == Tile::WALL)
                continue;

            ret.push_back(a);
        }

        // If knight moves are available, add them
        if (v.knight_cooldown == 0) {
            using enum LimpingKnightAction;
            for (LimpingKnightAction a: {
                     UP_UP_LEFT,
                     UP_UP_RIGHT,
                     DOWN_DOWN_LEFT,
                     DOWN_DOWN_RIGHT,
                     LEFT_LEFT_UP,
                     LEFT_LEFT_DOWN,
                     RIGHT_RIGHT_UP,
                     RIGHT_RIGHT_DOWN
                 }) {
                Position next_pos = perform_move(v.pos, a);
                if (next_pos.row >= _maze.size() ||
                    next_pos.col >= _maze[next_pos.row].size() ||
                    _maze[next_pos.row][next_pos.col] == Tile::WALL)
                    continue;
                ret.push_back(a);
            }
        }

        return ret;
    }

    State next_state(State v, const Action &a) const {
        int cooldown = v.knight_cooldown;

        // Check if the move is a knight move
        if (is_knight_move(a)) {
            cooldown = 7; // Reset cooldown after a knight move
        } else {
            // Decrease cooldown on adjacent moves or rest
            cooldown = (cooldown > 0) ? cooldown - 1 : 0;
        }

        // Return the updated state with new position and cooldown
        return {perform_move(v.pos, a), cooldown};
    }

    // Check if the current state is the target
    bool is_target(State v) const {
        return _maze[v.pos.row][v.pos.col] == Tile::EXIT;
    }

private:
    Maze _maze;
    Position _hero;
};
#endif


// A problem from lecture:
// You are given set of bottles of (possibly) different capacities,
// infinite supply of water and your task is to prepare exactly given
// amount of water in any of the bottles using as few operations
// as possible.
//
// The possible operations are:
// - Fill given bottle to full.
// - Empty given bottle.
// - Transfer water from bottle i to bottle j for i != j. Transfer stops
//   when i is empty or j is full whichever happens first.
#if TEST_BottleProblem
struct BottleProblem {
    struct State {
        std::vector<unsigned> water_in_bottles;

        // This will allow comparing states for using in sets or maps
        friend auto operator <=>(const State &, const State &) = default;
    };

    using Action = BottleOp;

    // Constructor to initialize target amount and bottle capacities
    BottleProblem(unsigned target_amount, const std::vector<unsigned> &bottle_capacities)
        : _target_amount(target_amount), _bottle_capacities(bottle_capacities) {
    }

    // Initial state is all bottles empty
    State initial_state() const {
        return {std::vector<unsigned>(_bottle_capacities.size(), 0)};
    }

    // Return all possible actions in the given state
    std::vector<Action> possible_actions(State v) const {
        std::vector<Action> actions;

        // Fill operations
        for (size_t i = 0; i < v.water_in_bottles.size(); ++i) {
            if (v.water_in_bottles[i] < _bottle_capacities[i]) {
                actions.push_back(Fill{.index = i});
            }
        }

        // Empty operations
        for (size_t i = 0; i < v.water_in_bottles.size(); ++i) {
            if (v.water_in_bottles[i] > 0) {
                actions.push_back(Empty{.index = i});
            }
        }

        // Transfer operations
        for (size_t i = 0; i < v.water_in_bottles.size(); ++i) {
            if (v.water_in_bottles[i] > 0) {
                // Only transfer from non-empty bottles
                for (size_t j = 0; j < v.water_in_bottles.size(); ++j) {
                    if (i != j && v.water_in_bottles[j] < _bottle_capacities[j]) {
                        actions.push_back(Transfer{.from = i, .to = j});
                    }
                }
            }
        }

        return actions;
    }

    // Apply action to get the next state
    State next_state(State v, const Action &a) const {
        std::visit(overloaded{
                       [&](const Fill &op) {
                           v.water_in_bottles[op.index] = _bottle_capacities[op.index];
                       },
                       [&](const Empty &op) {
                           v.water_in_bottles[op.index] = 0;
                       },
                       [&](const Transfer &op) {
                           unsigned transfer_amount = std::min(
                               v.water_in_bottles[op.from],
                               _bottle_capacities[op.to] - v.water_in_bottles[op.to]
                           );
                           v.water_in_bottles[op.from] -= transfer_amount;
                           v.water_in_bottles[op.to] += transfer_amount;
                       }
                   }, a);
        return v;
    }

    // Check if any bottle contains the target amount of water
    bool is_target(State v) const {
        return std::any_of(v.water_in_bottles.begin(), v.water_in_bottles.end(),
                           [&](unsigned amount) { return amount == _target_amount; });
    }

private:
    unsigned _target_amount;
    std::vector<unsigned> _bottle_capacities;
};
#endif


#ifndef __PROGTEST__

namespace maze_tests {

constexpr auto _ = Tile::EMPTY;
constexpr auto W = Tile::WALL;
constexpr auto X = Tile::EXIT;

struct Test {
  struct {
    size_t lost_hero, superstitious_hero, limping_knight;
  } solution;
  Position hero;
  Maze maze;
};

inline const Test TESTS[] = {
  { { 0, 0, 0 }, Position{0, 0}, Maze{ { X } } },
  { { 9, 13, 9 }, Position{ 0, 0 }, Maze{
    { _ },
    { _ },
    { _ },
    { _ },
    { _ },
    { _ },
    { _ },
    { _ },
    { _ },
    { X },
  } },
  { { 6, 6, 2 }, Position{ 1, 2 }, Maze{
    { W, _, _, _, W, },
    { W, _, _, W, X, },
    { _, _, _, W, _, },
    { _, _, _, _, _, },
    { _, W, _, W, _, },
  } },
  { { 2, 2, 2 }, Position{ 1, 2 }, Maze{
    { W, _, W, _, _, },
    { _, _, _, _, W, },
    { _, _, _, _, _, },
    { W, _, X, _, W, },
    { _, _, W, W, W, },
  } },
  { { 6, 6, 4 }, Position{ 3, 0 }, Maze{
    { W, _, _, X, _, },
    { _, _, _, _, _, },
    { W, _, _, _, W, },
    { _, _, W, W, _, },
    { _, _, _, W, W, },
  } },
  { { NO_SOLUTION, NO_SOLUTION, 1 }, Position{ 3, 2 }, Maze{
    { _, _, _, W, _, },
    { _, _, _, _, W, },
    { W, _, _, _, W, },
    { W, _, _, W, _, },
    { X, W, _, _, _, },
  } },
  { { 13, 19, 9 }, Position{ 4, 1 }, Maze{
    { W, W, _, W, _, _, _, _, _, _, _, W, W, _, _, },
    { _, _, _, _, W, W, W, _, _, W, W, _, _, _, _, },
    { _, _, _, _, _, W, _, W, _, _, W, _, _, _, _, },
    { W, _, _, W, _, W, _, _, _, _, W, W, W, _, W, },
    { _, _, _, _, W, W, _, _, _, _, _, _, W, _, _, },
    { _, _, W, _, W, _, _, W, _, _, _, _, _, _, _, },
    { _, _, _, W, _, _, W, W, W, _, _, _, W, _, _, },
    { _, W, _, _, _, _, _, _, _, W, _, _, _, _, W, },
    { W, _, W, _, _, W, _, W, _, _, _, W, _, _, _, },
    { _, W, _, _, _, _, _, _, _, X, W, _, W, _, W, },
  } },
  { { 6, 8, 4 }, Position{ 3, 5 }, Maze{
    { _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, },
    { _, _, _, _, W, W, _, _, _, _, _, _, _, W, _, },
    { X, _, _, _, _, _, _, W, _, _, _, _, _, _, _, },
    { W, _, _, _, W, _, _, _, _, _, W, _, _, W, _, },
    { _, W, _, _, _, _, W, _, _, _, _, W, W, _, W, },
    { W, _, _, _, _, _, _, W, W, _, _, _, _, W, _, },
    { _, _, _, _, W, _, W, W, W, _, _, _, W, _, _, },
    { W, _, W, _, _, W, _, _, W, W, _, W, W, W, _, },
    { _, _, _, _, W, _, W, W, W, _, _, W, _, _, _, },
    { _, _, _, _, _, _, _, _, _, _, W, _, _, _, _, },
  } },
  { { NO_SOLUTION, NO_SOLUTION, 3 }, Position{ 8, 6 }, Maze{
    { _, _, _, W, _, W, _, _, W, _, _, _, _, W, _, },
    { _, _, _, W, _, _, _, _, _, _, _, W, _, W, W, },
    { _, _, _, _, _, _, W, _, _, W, W, _, _, _, W, },
    { _, W, _, _, _, _, _, _, _, W, _, _, _, W, _, },
    { _, _, _, _, W, _, _, _, _, _, _, _, W, W, _, },
    { _, W, _, _, _, _, W, _, _, W, W, W, _, _, _, },
    { W, W, W, _, _, W, W, _, _, X, _, _, _, W, W, },
    { W, _, _, W, W, _, W, _, W, _, _, _, _, _, _, },
    { W, _, W, _, _, _, _, W, _, _, _, _, _, _, _, },
    { _, _, _, _, _, W, _, W, _, _, _, _, W, _, _, },
  } },
  { { 10, 12, 8 }, Position{ 7, 1 }, Maze{
    { W, _, _, _, _, _, _, _, _, W, W, W, _, _, _, },
    { W, _, _, W, _, W, _, _, W, _, _, _, _, _, _, },
    { W, W, _, _, _, _, _, _, W, _, W, _, W, _, W, },
    { _, _, _, _, _, _, _, _, _, W, _, W, W, _, W, },
    { _, W, W, _, _, _, _, _, X, _, _, _, _, W, _, },
    { W, _, _, _, _, W, W, _, _, _, _, _, _, _, W, },
    { _, W, _, _, W, _, _, _, W, _, _, _, _, W, _, },
    { _, _, _, W, _, _, _, _, _, _, _, W, _, _, W, },
    { _, _, _, _, _, W, _, _, W, W, _, W, _, W, _, },
    { W, _, _, W, _, _, W, W, W, W, _, _, _, _, _, },
  } },
  { { 46, 74, 36 }, Position{ 0, 27 }, Maze{
    { _, _, _, _, _, W, _, _, _, _, W, _, W, W, _, _, _, W, W, _, _, W, _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, _, W, _, _, _, _, W, _, _, W, _, _, _, },
    { _, _, _, _, W, W, W, _, _, _, _, _, W, _, _, _, _, _, _, _, W, _, _, W, W, _, W, _, _, _, W, _, _, _, _, _, W, _, _, _, _, _, _, _, W, _, _, W, _, _, },
    { _, W, _, _, _, _, W, W, _, _, _, _, W, _, W, W, W, _, _, _, _, W, _, W, W, _, _, _, _, _, _, _, W, W, W, _, _, W, _, W, _, _, _, W, W, _, _, W, _, _, },
    { _, W, _, W, _, _, W, _, _, _, _, _, W, W, _, _, _, _, _, _, _, W, _, _, W, W, _, _, _, _, W, W, _, W, _, _, W, _, _, _, W, _, W, _, _, _, _, _, W, _, },
    { _, W, W, _, _, _, _, _, _, _, _, W, _, _, W, W, _, _, _, _, W, _, W, _, W, _, _, _, W, _, _, _, _, _, W, _, _, W, _, _, _, _, _, _, W, _, _, _, _, _, },
    { _, _, _, W, W, _, _, _, _, _, _, _, W, W, W, W, _, W, _, _, _, _, W, W, _, _, _, _, _, _, _, _, _, W, _, _, _, W, _, _, _, W, W, _, W, _, W, W, _, _, },
    { W, _, _, _, _, _, _, _, _, _, _, W, W, _, _, W, _, W, W, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, _, _, W, _, W, _, W, _, W, W, W, _, },
    { W, W, _, _, W, _, _, W, _, _, W, _, _, W, W, _, W, _, _, _, _, _, _, W, W, _, _, _, W, _, W, _, _, _, _, W, _, _, _, _, _, W, _, _, _, W, _, W, _, _, },
    { _, W, W, W, _, _, W, W, _, _, W, W, _, _, _, _, _, _, _, _, _, W, W, _, _, W, _, _, W, _, _, _, _, W, _, _, _, _, W, _, _, W, _, _, W, _, W, _, _, _, },
    { _, _, _, _, W, _, _, W, W, W, _, W, _, _, _, W, _, W, _, _, _, W, _, W, _, _, _, W, W, W, _, W, _, _, _, _, W, _, _, _, _, W, _, _, _, _, _, _, W, W, },
    { _, _, _, _, W, _, W, _, _, _, W, _, _, _, _, _, W, _, W, W, W, _, W, W, _, W, W, W, W, W, _, W, _, _, _, _, _, W, W, W, _, _, _, _, _, _, W, _, _, _, },
    { W, _, _, _, _, W, _, _, W, _, _, _, _, _, W, _, W, W, _, _, _, _, _, W, _, _, _, _, _, W, W, _, _, W, W, _, _, _, _, _, _, _, _, _, _, _, _, W, _, _, },
    { _, _, _, _, W, _, W, _, _, W, _, _, _, _, _, W, _, W, _, W, _, _, _, W, _, _, _, _, _, _, _, _, _, W, W, _, _, W, _, _, _, _, _, _, _, W, _, _, W, _, },
    { _, _, _, _, _, _, W, W, W, W, _, W, _, _, _, W, _, _, _, _, W, W, _, _, _, _, _, _, _, W, _, W, _, _, _, W, W, W, _, _, _, _, _, _, _, _, _, _, _, _, },
    { _, _, W, _, _, W, _, _, _, _, _, _, _, _, _, W, _, _, W, W, _, W, W, _, _, _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, _, _, W, _, _, _, _, W, W, W, },
    { _, _, _, W, W, _, _, _, _, _, _, _, _, W, W, _, _, _, _, _, _, _, _, _, _, _, W, _, _, _, W, W, W, _, _, W, _, _, W, _, _, _, _, _, W, _, _, _, W, W, },
    { W, _, _, W, _, _, _, W, _, W, _, _, _, _, W, _, _, _, _, _, _, _, W, _, _, _, W, W, W, _, W, _, W, _, _, W, W, _, W, _, W, _, _, _, _, _, W, _, _, _, },
    { _, W, W, _, _, W, _, _, _, _, _, _, W, _, _, W, _, _, W, W, _, _, _, _, _, W, W, _, _, W, _, W, _, _, _, _, W, _, _, _, W, W, W, _, _, _, _, _, _, W, },
    { W, W, _, _, W, _, _, W, _, W, _, W, _, _, _, _, _, W, W, _, _, _, _, W, _, _, _, _, W, _, _, _, _, _, _, _, _, W, _, _, _, _, _, W, _, _, _, _, _, _, },
    { _, _, _, _, _, _, W, W, W, _, _, _, _, _, W, _, _, W, _, _, W, _, _, _, _, _, _, W, W, W, _, _, _, W, W, W, _, _, _, _, _, _, W, W, _, _, _, _, _, _, },
    { W, _, W, _, _, W, W, W, W, W, _, _, W, _, _, _, _, _, _, _, _, W, _, _, W, _, _, _, _, _, W, _, W, W, _, _, _, _, W, _, _, _, W, _, _, _, W, _, _, _, },
    { _, _, _, _, W, _, _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, _, W, W, _, _, _, _, W, _, _, W, _, _, _, W, W, _, _, _, _, _, _, W, W, _, _, },
    { _, _, _, _, W, _, W, W, _, _, _, _, W, _, _, _, _, W, _, _, _, _, _, _, W, _, _, _, _, W, W, _, _, _, _, W, _, W, _, W, W, W, _, _, _, W, _, _, W, _, },
    { _, W, _, _, _, _, _, W, _, _, _, _, W, _, _, _, W, _, _, _, _, W, _, _, _, _, _, _, _, W, _, _, W, _, _, W, _, _, _, _, _, _, _, _, _, _, _, _, W, _, },
    { _, _, _, _, _, _, _, W, _, _, _, _, _, _, W, _, _, _, W, _, _, _, _, W, _, _, _, W, _, _, _, _, _, _, _, _, W, _, W, _, _, _, _, _, W, W, _, _, W, _, },
    { _, W, _, W, W, _, W, W, W, _, W, _, _, _, W, _, _, _, _, _, _, _, _, _, W, _, _, _, W, _, W, W, W, _, _, _, _, _, _, _, _, _, _, _, W, W, _, W, _, _, },
    { _, W, _, _, W, _, _, _, _, _, _, _, W, W, _, _, W, _, _, _, W, _, _, _, _, W, _, _, _, _, _, W, W, _, _, _, W, _, _, W, _, W, W, _, W, W, _, X, _, _, },
    { W, _, _, W, _, W, _, W, _, _, _, W, _, W, _, _, _, _, _, _, _, _, W, _, _, W, _, _, W, _, W, W, _, _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, W, W, },
    { _, _, _, _, _, W, _, _, _, _, W, _, _, _, W, _, W, _, _, W, _, _, _, _, _, _, W, _, _, _, _, _, W, _, _, _, _, _, _, W, _, W, W, _, _, _, _, _, W, W, },
    { _, W, _, W, _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, W, W, _, _, W, _, _, _, W, _, W, W, _, _, W, W, _, _, W, W, _, _, W, _, W, },
  } },
  { { NO_SOLUTION, NO_SOLUTION, 12 }, Position{ 17, 34 }, Maze{
    { _, _, _, _, _, _, _, _, W, _, W, W, _, _, W, _, _, _, _, _, _, _, W, W, _, W, _, _, _, _, _, _, _, W, _, W, W, W, _, _, _, W, _, W, _, _, _, _, _, _, },
    { _, _, _, _, _, W, _, _, _, _, W, _, _, _, _, _, W, _, W, W, _, _, _, _, _, _, _, _, W, W, W, _, _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, W, W, },
    { W, W, W, _, _, _, _, _, W, _, W, _, W, _, _, _, _, W, _, _, W, _, W, W, _, _, W, _, W, _, _, W, W, _, _, W, _, _, _, _, _, W, _, _, W, _, _, _, _, W, },
    { _, _, W, _, W, W, W, _, W, W, W, _, _, W, _, _, _, _, _, _, _, W, W, W, _, W, W, _, W, _, _, _, _, _, _, _, W, _, W, _, _, _, _, _, W, _, _, _, _, W, },
    { _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, W, _, _, W, _, _, W, _, _, _, W, W, W, _, W, _, W, _, W, _, _, _, _, _, _, W, _, },
    { _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, _, _, _, _, W, _, _, _, W, W, _, _, _, W, _, W, W, W, W, _, W, _, W, _, _, _, _, _, _, W, _, },
    { W, W, _, W, _, _, _, _, _, _, _, _, _, W, _, _, W, _, _, W, _, _, _, W, _, _, _, _, W, _, _, _, W, W, _, W, _, _, W, W, _, _, W, _, W, W, _, _, _, _, },
    { _, _, W, W, W, W, W, _, _, _, _, _, W, _, W, W, _, _, W, W, _, _, _, _, W, W, _, _, W, W, _, _, _, _, W, W, _, _, _, W, _, _, _, _, W, _, W, _, _, W, },
    { W, _, W, _, _, _, _, _, W, _, _, _, W, _, _, W, _, _, _, _, _, _, _, W, _, _, _, _, _, _, W, _, _, _, W, _, _, _, _, _, W, _, W, W, W, _, _, _, _, _, },
    { _, _, W, _, _, _, W, _, _, W, _, _, _, _, W, W, _, _, _, _, _, _, _, _, _, _, W, _, _, W, _, _, _, W, W, _, _, _, _, W, W, W, _, _, _, W, W, _, _, _, },
    { _, _, W, _, W, W, W, W, W, _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, W, _, _, W, W, W, W, _, _, _, _, W, W, W, _, _, W, _, },
    { _, _, W, _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, W, _, W, W, _, _, _, W, _, _, W, W, _, W, _, _, _, W, _, _, W, _, _, _, _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, _, _, W, W, _, W, W, _, W, W, _, W, W, W, _, _, W, _, W, W, _, W, W, _, _, _, _, _, _, _, W, _, _, _, W, _, _, W, _, W, _, _, _, _, },
    { _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, _, _, W, _, _, _, _, _, _, W, _, _, _, _, W, _, W, _, _, _, _, _, _, W, W, _, W, W, _, },
    { _, W, _, _, _, _, _, W, W, W, _, _, _, W, W, _, _, _, _, _, W, _, _, W, W, _, _, _, _, W, _, _, _, _, _, _, _, W, _, _, _, W, W, _, _, _, W, _, W, _, },
    { _, _, _, _, _, W, _, _, _, W, _, W, _, _, W, W, W, _, _, _, W, _, _, _, _, _, W, _, _, W, _, _, _, W, _, W, _, _, _, W, _, _, W, _, _, _, _, _, _, _, },
    { _, W, _, W, W, _, _, W, _, _, _, W, _, _, _, _, W, _, _, W, _, _, _, W, _, _, W, _, _, W, _, W, _, _, _, W, W, _, _, _, _, W, _, _, W, _, W, _, _, _, },
    { _, _, _, W, _, _, W, _, _, _, _, _, _, W, _, W, _, W, W, _, W, _, _, W, _, _, _, _, _, _, _, _, _, _, _, W, _, _, W, W, _, _, _, _, W, _, _, W, _, _, },
    { _, _, _, W, _, _, _, W, _, _, _, _, W, _, _, _, _, W, W, _, _, _, _, _, W, _, W, _, _, _, W, W, _, _, _, _, W, W, _, _, _, W, _, _, _, _, _, _, _, _, },
    { W, W, _, W, W, W, _, _, _, _, W, W, _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, W, _, W, _, _, _, W, _, W, W, _, _, _, W, W, _, _, _, _, _, W, },
    { _, _, _, _, _, _, _, W, _, _, _, W, _, _, W, _, _, _, _, _, _, _, W, W, W, _, W, _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, W, W, _, W, _, _, _, _, },
    { W, W, _, W, _, _, _, _, _, _, _, W, _, _, _, W, _, _, _, _, W, _, W, _, _, W, _, _, W, _, W, W, _, _, _, _, _, W, _, W, W, W, W, _, _, _, _, _, W, _, },
    { _, _, _, _, W, W, W, _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, W, W, _, _, _, W, _, _, },
    { _, _, _, W, _, _, _, _, _, W, _, _, _, _, W, _, _, _, _, _, W, _, _, _, _, _, _, W, _, _, _, W, W, _, W, _, W, W, _, _, W, _, _, W, _, _, _, _, _, _, },
    { _, W, _, _, _, W, _, _, _, _, _, _, _, _, W, _, _, W, W, _, _, _, _, W, W, _, _, W, _, _, _, _, _, W, _, _, _, _, _, _, _, _, W, _, W, _, _, _, _, _, },
    { _, _, _, W, _, _, _, _, W, W, _, _, W, W, _, W, _, _, _, _, _, W, _, _, _, _, _, W, W, _, _, _, W, _, W, W, _, _, _, _, _, W, _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, _, W, W, _, _, _, _, W, _, _, _, _, W, _, _, W, _, W, W, _, _, _, _, W, _, W, _, _, W, _, _, W, _, _, W, _, W, _, _, W, _, W, _, W, },
    { _, W, _, _, W, W, _, W, _, W, _, _, _, _, _, _, W, W, _, _, _, _, _, _, _, _, _, _, W, _, _, _, W, _, _, W, _, _, _, _, _, W, _, W, _, _, _, _, W, _, },
    { W, _, _, _, _, _, _, _, W, _, _, _, _, _, _, W, _, _, _, _, W, W, _, _, _, _, _, _, _, _, _, _, W, _, _, _, W, W, W, _, _, _, _, W, _, _, _, _, _, W, },
    { _, _, _, _, _, _, W, _, W, _, _, _, _, _, W, _, _, _, _, W, _, _, W, _, W, W, _, W, W, _, _, _, _, W, _, W, X, _, W, W, _, _, W, _, W, _, _, W, _, W, },
  } },
  { { 19, 37, 9 }, Position{ 14, 29 }, Maze{
    { _, W, _, _, W, _, _, _, _, _, W, _, _, W, W, W, _, _, _, W, _, _, _, _, _, _, _, _, _, W, _, _, _, _, W, _, W, W, _, _, W, W, W, _, W, W, W, _, _, _, },
    { _, _, _, _, _, W, _, _, _, W, W, W, W, W, _, _, _, _, _, _, _, _, _, _, W, _, _, _, _, _, W, W, _, _, _, _, _, W, W, W, _, W, _, _, W, _, _, _, _, W, },
    { W, W, _, _, W, _, _, W, _, _, _, W, _, _, _, _, _, _, W, _, _, W, _, _, _, W, _, _, _, _, _, _, W, _, _, W, _, _, _, _, _, _, _, W, _, _, _, _, _, _, },
    { W, _, _, _, W, _, _, _, _, _, _, _, W, _, _, _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, W, _, _, _, W, W, _, _, _, W, _, W, _, _, _, _, W, _, _, },
    { _, _, W, _, W, _, W, _, _, _, _, W, _, _, _, W, _, _, _, W, _, _, W, _, _, _, W, _, W, _, _, W, _, _, _, W, _, _, _, _, _, _, W, _, _, _, _, _, _, _, },
    { _, W, _, _, _, _, W, W, _, _, _, _, _, _, _, _, _, _, _, W, W, _, _, _, _, _, W, _, _, W, _, _, W, W, _, W, _, _, _, _, _, _, W, _, _, _, _, _, _, W, },
    { _, W, _, W, W, _, W, _, _, _, W, _, _, _, _, W, _, W, W, _, W, W, W, _, _, W, _, W, W, _, _, W, W, _, W, _, _, _, _, W, _, _, W, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, _, _, W, _, _, W, _, _, _, _, W, _, _, _, _, W, _, _, _, _, W, _, _, _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, W, _, _, W, _, },
    { W, _, _, _, _, W, _, W, _, _, _, _, _, _, _, W, _, _, _, _, W, W, _, W, _, W, _, _, W, _, W, _, _, _, _, _, _, _, W, W, W, _, _, W, W, _, W, _, _, _, },
    { _, _, _, _, _, _, _, W, _, _, _, _, W, _, _, _, _, W, _, _, W, W, _, _, _, _, _, W, _, _, _, _, _, W, _, _, W, _, _, W, _, _, _, _, W, W, _, _, W, _, },
    { W, _, W, _, W, _, W, W, _, _, _, _, _, _, _, W, _, _, W, _, W, W, W, _, _, _, _, _, _, W, _, _, _, _, W, W, _, _, _, _, _, W, _, _, _, _, _, _, _, _, },
    { _, W, W, _, _, _, W, _, _, _, _, W, W, W, W, _, _, _, _, _, _, W, W, _, _, W, _, _, _, _, _, W, W, _, W, W, _, W, W, W, _, _, _, _, _, _, _, _, _, W, },
    { W, W, W, _, _, _, _, _, _, W, W, _, _, W, _, _, _, W, W, _, _, _, W, _, W, _, W, _, _, _, _, W, _, _, _, W, _, _, _, W, _, _, _, _, _, _, _, W, _, _, },
    { _, W, _, _, W, W, _, _, _, W, W, W, _, _, _, _, _, _, _, _, W, W, W, _, _, _, _, _, _, W, _, _, _, _, _, _, W, W, _, _, _, W, W, W, _, _, W, _, _, W, },
    { _, _, _, W, _, W, _, W, _, _, _, W, _, W, W, W, _, _, W, _, W, _, W, _, _, _, _, W, _, _, _, _, W, W, _, _, _, W, _, W, _, _, _, _, W, _, _, _, W, _, },
    { _, W, _, _, _, _, _, _, _, _, W, W, _, W, _, _, W, _, _, _, W, _, W, W, _, _, W, _, W, W, W, W, W, _, W, _, _, _, _, W, _, _, _, W, _, _, _, _, _, _, },
    { W, _, _, W, _, _, _, _, _, _, _, W, W, _, W, _, _, W, W, _, _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, W, W, W, W, _, _, _, W, _, _, _, _, _, _, },
    { _, _, _, W, _, _, _, W, _, W, W, _, W, W, _, _, W, _, W, W, _, _, _, _, W, _, _, _, _, _, _, _, _, W, _, W, _, _, _, _, _, _, W, _, _, W, W, _, _, W, },
    { _, _, _, _, _, _, _, _, _, W, _, W, _, _, _, W, W, _, _, _, W, _, _, _, W, _, _, _, _, _, W, _, _, _, W, _, _, _, _, _, _, _, _, _, _, _, _, W, _, _, },
    { _, _, _, _, W, _, W, _, W, _, W, W, _, _, _, W, _, _, _, _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, W, _, W, _, W, _, _, _, _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, _, _, W, W, _, W, _, _, _, _, _, _, _, W, _, _, _, _, W, W, W, _, W, _, _, W, _, _, X, _, _, W, _, _, _, _, _, W, _, _, _, _, _, W, },
    { W, _, _, _, _, _, _, _, _, _, W, _, _, W, _, W, _, _, _, W, _, _, _, _, _, _, _, _, _, W, _, _, _, W, W, W, _, _, _, _, W, W, _, _, _, _, W, _, _, _, },
    { _, _, _, _, _, _, _, _, W, _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, W, _, _, _, _, W, W, W, _, W, _, W, _, _, _, _, W, _, _, W, _, _, _, W, _, _, },
    { _, _, _, W, _, _, _, W, _, _, _, _, _, _, _, W, _, W, W, _, _, W, _, _, W, W, _, _, _, _, W, W, _, _, _, _, W, W, _, _, _, _, _, _, _, _, _, W, _, _, },
    { _, W, W, _, _, _, _, _, W, W, W, _, W, W, _, W, W, _, _, _, _, _, _, _, _, _, _, _, W, _, _, _, _, _, W, W, _, _, _, W, _, _, W, _, W, _, _, _, _, _, },
    { W, _, W, _, _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, _, _, W, _, W, _, _, _, _, _, _, W, W, W, W, _, _, _, W, W, },
    { _, W, W, _, _, W, _, _, _, _, _, _, _, _, W, W, _, W, _, _, _, W, W, W, _, _, _, _, _, _, W, _, _, _, W, _, W, _, _, _, _, _, _, W, W, W, W, _, _, W, },
    { _, _, _, W, _, _, W, W, _, _, _, _, _, _, _, _, _, _, _, _, _, W, _, _, W, _, _, _, _, W, W, _, _, W, _, _, _, _, W, _, _, W, _, _, W, W, W, _, _, W, },
    { _, _, W, _, _, _, W, W, _, _, W, _, _, W, W, _, _, _, _, _, W, _, _, _, _, W, _, W, _, _, _, _, _, _, W, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, },
    { _, _, _, W, _, _, _, _, _, _, _, _, W, _, _, W, _, W, _, W, W, _, _, _, _, _, _, _, _, _, _, _, _, W, _, W, _, _, _, W, _, _, _, _, _, _, _, _, _, W, },
  } },
  { { 23, 27, 17 }, Position{ 22, 4 }, Maze{
    { _, W, W, W, _, _, _, _, W, _, _, _, _, _, _, W, _, W, W, _, W, _, _, W, _, _, W, W, _, W, W, _, W, _, W, W, W, _, _, _, W, _, W, _, _, _, W, _, _, _, },
    { _, _, _, _, _, _, _, W, W, _, _, W, _, _, W, _, W, _, W, _, _, _, _, _, _, W, W, _, W, W, W, _, _, W, _, W, W, _, W, W, _, _, W, _, _, W, W, _, W, _, },
    { _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, _, _, _, _, W, _, W, _, _, _, _, W, _, _, _, _, _, _, _, _, W, _, W, _, _, W, _, _, _, W, },
    { W, _, _, _, _, W, _, W, _, _, _, _, W, W, W, _, _, W, _, _, _, _, _, _, W, _, _, _, _, W, _, _, W, _, W, W, W, _, _, _, _, _, _, _, W, W, _, _, _, W, },
    { _, _, _, _, _, _, _, _, _, _, _, W, W, _, _, _, _, W, _, W, _, _, _, _, _, _, _, _, W, _, _, _, _, _, W, _, W, _, _, _, _, _, _, W, _, _, _, _, W, _, },
    { W, _, _, _, _, _, _, W, W, W, _, _, _, _, _, _, W, _, _, _, _, W, W, _, W, _, _, _, _, W, _, W, _, _, _, _, W, W, _, W, _, W, _, _, _, _, _, _, _, _, },
    { W, W, _, W, _, _, _, W, _, _, W, _, _, _, _, W, _, _, _, _, _, W, _, _, _, W, _, _, _, _, _, _, _, _, W, _, W, _, _, _, _, W, W, W, _, _, W, _, _, W, },
    { _, _, W, _, W, W, W, _, W, _, _, _, _, _, W, _, W, _, _, _, _, _, W, _, _, _, _, _, _, W, _, W, _, W, W, _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, _, W, _, W, W, _, _, W, _, _, _, _, W, _, W, _, _, _, _, _, _, _, _, _, _, _, W, W, _, _, _, _, _, _, _, },
    { _, _, _, _, _, W, _, _, W, W, W, W, _, W, _, _, _, W, _, _, W, W, W, _, _, _, _, _, W, W, _, _, W, _, _, W, _, _, _, W, W, _, W, _, _, _, _, _, W, W, },
    { _, _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, W, _, W, _, W, W, _, _, _, _, W, W, _, _, _, _, _, _, _, _, W, W, W, _, W, _, _, _, _, _, W, _, _, _, },
    { _, W, _, W, _, _, W, _, _, _, _, _, _, _, X, _, W, _, _, _, _, W, W, W, _, _, _, W, _, W, _, _, _, _, _, W, _, _, _, _, _, _, W, W, W, _, _, _, W, W, },
    { _, _, _, _, _, _, _, _, _, _, W, _, _, _, W, _, _, _, _, _, _, _, _, W, W, _, W, W, _, _, _, W, W, W, W, W, _, _, _, _, _, _, W, _, _, _, _, W, _, _, },
    { _, _, _, _, _, _, W, _, _, _, _, W, W, W, W, _, _, _, _, _, _, W, _, W, W, _, _, _, _, _, _, W, _, _, _, _, _, W, _, _, _, _, _, _, _, _, W, _, _, _, },
    { _, _, _, _, _, _, W, _, _, W, W, _, _, _, _, W, _, _, _, _, W, _, _, W, _, W, _, _, _, _, _, W, _, W, _, _, W, _, W, _, _, W, W, W, _, _, _, _, W, _, },
    { _, _, W, _, _, _, _, _, _, _, _, W, _, _, W, W, _, W, _, _, _, _, W, _, W, _, _, _, _, _, W, W, _, _, _, W, W, W, _, W, _, W, W, _, W, _, _, _, _, W, },
    { _, _, _, _, W, W, W, W, W, _, _, _, _, _, W, W, W, W, W, _, _, _, _, _, W, _, _, _, W, W, _, _, _, W, _, W, W, _, W, _, _, _, _, _, _, _, _, _, _, _, },
    { _, _, _, _, W, _, W, _, W, _, _, _, _, _, W, _, W, _, _, W, _, _, _, _, W, _, W, _, _, W, _, _, _, _, W, _, W, _, W, _, _, _, W, _, _, _, W, _, _, _, },
    { _, _, _, _, _, _, _, _, _, _, _, _, W, _, W, _, _, _, _, _, W, _, _, _, _, _, _, W, W, W, _, _, W, _, _, _, _, _, _, W, _, _, W, W, W, W, _, _, W, _, },
    { W, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, W, _, W, W, _, _, _, _, _, W, W, _, _, _, },
    { _, _, W, W, W, W, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, W, W, _, _, W, _, W, _, _, W, W, W, _, _, _, _, _, _, _, _, _, _, _, W, W, W, _, _, },
    { W, W, _, _, _, _, _, _, W, W, _, _, W, _, _, W, _, _, W, _, _, _, _, _, _, _, _, W, _, _, _, _, W, W, _, _, _, _, W, _, _, _, _, _, _, _, _, _, _, W, },
    { _, _, W, W, _, W, _, W, _, _, _, W, _, _, _, _, _, _, _, W, _, _, W, _, _, _, _, _, W, W, W, _, _, _, _, _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, },
    { W, _, _, W, _, _, _, W, _, W, _, _, W, _, _, _, _, _, _, W, W, _, _, _, _, _, _, _, W, _, _, W, W, _, _, _, W, _, _, _, _, _, _, _, _, _, _, _, _, W, },
    { W, _, _, _, _, W, W, W, W, _, _, _, _, _, W, _, W, _, _, _, _, _, _, W, W, _, W, _, W, W, W, _, _, _, _, W, W, _, W, _, W, _, _, W, W, _, _, _, W, _, },
    { _, _, _, _, _, W, _, _, _, _, _, _, _, _, _, W, _, _, _, _, W, W, _, _, _, _, _, W, _, W, _, _, W, W, W, _, _, _, _, _, W, _, _, _, _, _, W, _, _, _, },
    { W, _, _, _, W, W, _, W, _, _, _, _, _, _, _, _, _, _, _, W, _, _, _, _, _, W, _, _, W, W, _, W, _, W, _, W, _, _, _, W, _, _, _, _, W, _, _, _, _, _, },
    { W, W, _, W, W, _, _, _, _, _, W, W, _, _, _, _, W, _, _, _, W, _, _, _, W, W, W, _, _, _, W, _, W, W, _, W, _, _, _, _, W, W, _, _, _, _, _, _, _, _, },
    { _, W, _, _, _, _, _, _, _, W, W, _, _, _, _, W, _, _, _, W, W, W, _, W, _, _, _, W, _, W, W, W, W, _, W, _, W, _, _, _, _, W, W, _, _, W, W, W, _, _, },
    { W, _, W, _, W, _, _, _, W, W, _, W, _, _, W, W, _, W, _, _, _, _, _, _, _, _, _, _, W, _, _, W, _, W, _, _, _, _, _, W, _, W, _, _, _, W, _, W, _, _, },
  } },
};

} // namespace maze_tests


namespace bottle_tests {

struct Test {
  size_t solution;
  unsigned target_amount;
  std::vector<unsigned> bottles;
};

inline const Test TESTS[] = {
  { 800, 600, { 1000, 1 } },
  { 126, 12, { 78, 61, } },
  { NO_SOLUTION, 53, { 7, 39, } },
  { 24, 5, { 82, 71, } },
  { NO_SOLUTION, 71, { 20, 29, } },
  { NO_SOLUTION, 55, { 27, 93, } },
  { 4, 29, { 75, 70, 24, } },
  { 14, 49, { 33, 13, 94, } },
  { NO_SOLUTION, 64, { 18, 56, 49, } },
  { NO_SOLUTION, 7, { 57, 15, 54, } },
  { 1, 7, { 36, 22, 4, 7, 26, } },
  { 6, 3, { 15, 40, 71, 76, 4, } },
  { 2, 29, { 89, 78, 63, 12, 99, 52, 41, } },
  { 4, 78, { 50, 20, 28, 57, 93, 9, 11, } },
  { 4, 31, { 21, 38, 77, 100, 1, 20, 56, 11, } },
  { 3, 34, { 74, 6, 60, 23, 25, 15, 60, 45, } },
  { 4, 21, { 49, 69, 57, 12, 93, 43, 77, 53, 46, 18, } },
};

}; // namespace bottle_tests

template < bool move_up, bool knight >
void check_hero(size_t ref, const Maze& maze, Position hero,
  const std::optional<std::vector<std::conditional_t<knight, LimpingKnightAction, Direction>>>& sol_
) {
  if (ref == NO_SOLUTION) {
    CHECK(!sol_, "Got solution but none should exist\n");
    return;
  }

  CHECK(sol_, "Solution exists but got empty optional\n");
  const auto &sol = sol_.value();

  Tile tile = maze[hero.row][hero.col];

  size_t move_cnt = 1;
  size_t jump_cd = 0;
  for (auto d : sol) {
    hero = perform_move(hero, d);

    if constexpr(move_up) {
      if (move_up && move_cnt++ % 5 == 0)
        CHECK(d == Direction::UP, "Every 5th move must be up!\n");
    }

    if constexpr(knight) {
      if (jump_cd > 0) jump_cd--;
      CHECK(jump_cd == 0 || !is_knight_move(d), "Hero cannot move as knight yet\n");
      if (is_knight_move(d)) jump_cd = 8;
    }

    CHECK(hero.row < maze.size(), "Row out of range\n");
    const auto& row = maze[hero.row];
    CHECK(hero.col < row.size(), "Col out of range\n");
    tile = row[hero.col];

    CHECK(tile == Tile::EMPTY || tile == Tile::EXIT, "Cannot step on wall\n");
  }

  CHECK(tile == Tile::EXIT, "Not at exit\n");
  CHECK(sol.size() == ref, "Expected %zu moves but got solution with %zu moves\n",
    ref, sol.size());
}


void check_bottles(size_t ref, unsigned target, const std::vector<unsigned>& bottles,
  const std::optional<std::vector<BottleOp>>& sol_
) {
  if (ref == NO_SOLUTION) {
    CHECK(!sol_, "Got solution but none should exist\n");
    return;
  }

  CHECK(sol_, "Solution exists but got empty optional\n");
  const auto &sol = sol_.value();

  std::vector<unsigned> water(bottles.size(), 0);

  for (auto op : sol) std::visit(overloaded{
    [&](Fill op) {
      CHECK(op.index < bottles.size(), "Bottle out of range\n");
      water[op.index] = bottles[op.index];
    },
    [&](Empty op) {
      CHECK(op.index < bottles.size(), "Bottle out of range\n");
      water[op.index] = 0;
    },
    [&](Transfer op) {
      CHECK(op.from < bottles.size(), "Source bottle out of range\n");
      CHECK(op.to < bottles.size(), "Target bottle out of range\n");
      CHECK(op.to != op.from, "Source and target bottles are the same\n");
      unsigned a = std::min(bottles[op.to] - water[op.to], water[op.from]);
      water[op.to] += a;
      water[op.from] -= a;
    }
  }, op);

  bool found_bottle = false;
  for (unsigned b : water) if (b == target) found_bottle = true;
  CHECK(found_bottle, "No bottle with target amount of water\n");

  CHECK(sol.size() == ref, "Expected %zu moves but got solution with %zu moves\n",
    ref, sol.size());
}

int main() {
# define T(x) (std::cout << #x << (x ? ": enabled" : ": disabled") << std::endl)
  T(TEST_SuperstitiousHero);
  T(TEST_LimpingKnight);
  T(TEST_BottleProblem);
# undef T

  size_t i = 0;
  for (auto&& t : maze_tests::TESTS) {
    i++;
    std::string where;
    try {
      {
        where = std::to_string(i) + " lost_hero";
        auto sol = solve(LostHero(t.maze, t.hero));
        check_hero<false, false>(t.solution.lost_hero, t.maze, t.hero, sol);
      }

#     if TEST_SuperstitiousHero
      {
        where = std::to_string(i) + " superstitious_hero";
        auto sol = solve(SuperstitiousHero(t.maze, t.hero));
        check_hero<true, false>(t.solution.superstitious_hero, t.maze, t.hero, sol);
      }
#     endif

#     if TEST_LimpingKnight
      {
        where = std::to_string(i) + " limping_knight";
        auto sol = solve(LimpingKnight(t.maze, t.hero));
        check_hero<false, true>(t.solution.limping_knight, t.maze, t.hero, sol);
      }
#     endif
    } catch (const TestFailed& e) {
      std::cout << "Test " << where << " failed: " << e.what() << std::endl;
      return 1;
    }
  }

  size_t j = 0;
# if TEST_BottleProblem
  for (auto&& t : bottle_tests::TESTS) {
    j++;
    try {
      auto sol = solve(BottleProblem(t.target_amount, t.bottles));
      check_bottles(t.solution, t.target_amount, t.bottles, sol);
    } catch (const TestFailed& e) {
      std::cout << "Test " << j << " failed: " << e.what() << std::endl;
      return 1;
    }
  }
# endif

  std::cout << "All " << i + j << " tests passed." << std::endl;
}

#endif


