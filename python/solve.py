"""Solves an instance.

Modify this file to implement your own solvers.

For usage, run `python3 solve.py --help`.
"""

import argparse
import math
from pathlib import Path
from typing import Callable, Dict

from instance import Instance
from point import Point # Added this in imports
from solution import Solution
from file_wrappers import StdinFileWrapper, StdoutFileWrapper


def solve_naive(instance: Instance) -> Solution:
    return Solution(
        instance=instance,
        towers=instance.cities,
    )

def algo(instance: Instance) -> Solution:
    print("in algo")
    D = instance.grid_side_length
    N = len(instance.cities)
    R_s = instance.coverage_radius
    R_p = instance.penalty_radius
    print("D", D)
    print("N", N)
    print("R_s", R_s)
    print("R_p", R_p)
    print("===================")

    # create P (potential tower) boolean arrays
    potential_towers = [[([False for _ in range(N)], Point(j, k)) for j in range(D)] for k in range(D)]
    for x in range (D):
        for y in range(D):
            for c in range(N):
                city_dist_from_tower = instance.cities[c].distance_sq(potential_towers[x][y][1])
                if city_dist_from_tower <= R_s ** 2:
                    potential_towers[x][y][0][c] = True
    
    # pick Pi with the most Ts
    cities_covered = [[(potential_towers[x][y][0], sum(potential_towers[x][y][0])) for x in range(D)] for y in range(D)]
    max_cities_covered = 0
    for x in range (D):
        for y in range (D):
            max_cities_covered = max(max_cities_covered, cities_covered[x][y])
            
    # holding all Pis with most cities covered
    potential_start_towers = set()
    for x in range (D):
        for y in range(D):
            if cities_covered[x][y][1] == max_cities_covered:
                potential_start_towers.add(Point(x, y), cities_covered[x][y][0])

    # are all cities visited?
    def checker(visited_cities, tower):
        for is_visited in range(N):
            visited_cities[is_visited] = tower[is_visited][1] or visited_cities[is_visited]

        return all(visited_cities), visited_cities

    # return a set of Pis that have least overlap with input Pi
    def least_overlap(tower, ans):
        ret_least_overlaps = set()
        min_overlap = float("inf")
        for x in range(D):
            for y in range(D):
                count = 0
                exists = True
                for t in ans:
                    exists = exists and (x == t[0][0] and y == t[0][1])
                if not exists:
                #if not (x == tower[0][0] and y == tower[0][1]):
                    for c in range(N):
                        if tower[1][c] and potential_towers[x][y][0][c]:
                            count += 1
                    min_overlap = min(min_overlap, count)

        for x in range(D):
            for y in range(D):
                if not (x == tower[0][0] and y == tower[0][1]) and sum(potential_towers[x][y][0]) == min_overlap:
                    ret_least_overlaps.add(potential_towers[x][y])
        return ret_least_overlaps

    # greedy 
    possible_solutions = []
    for tower in potential_start_towers:
        # visited cities
        visited_cities = [False for _ in range(N)]
        # set of towers in this iteration
        ans = set()
        ans.add(tower)
        #noting that some cities are visited with possible start tower
        check, visited_cities = checker(visited_cities, tower)
        if check:
           continue 
        
        most_recent_tower = tower
        while not check:
            next_towers = least_overlap(tower, ans)
            max_dist = 0
            next_tower = 0; # next tower that is furthest from tower
            if len(next_towers) > 1:
                for tow in next_towers:
                    temp = max_dist
                    max_dist = max(max_dist, tow[0].distance_sqr(most_recent_tower[0]))
                    if not temp == max_dist:
                        next_tower = tow 
            else:
                next_tower = next_towers.pop()
            ans.add(next_tower)
            most_recent_tower = next_tower
            check, visited_cities = checker(tower, next_tower)
        possible_solutions.append(ans)


    #check penalties
    def penalty(M): #M has 
        #to implement
        def count_overlap(tower):
            count = 0
            for i in M:
                #checking if i is itself
                dist = i[0].distance_sqr(tower[0])
                if dist != 0 and dist < R_p:
                    count += 1
            return count        
        return_val = 0
        for tower in M:
            w_j = count_overlap(tower)
            return_val += (170 * math.e ** (0.17 * w_j))
        return return_val
    min_penalty = float("inf")
    solution = 0
    for ans in possible_solutions:
        temp = min_penalty
        min_penalty = min(min_penalty, penalty(ans))
        if not temp == min_penalty:
            solution = ans
        '''for tower in ans:
            temp = min_penalty
            min_penalty = min(min_penalty, penalty(ans))
            if not temp == min_penalty:
                solution = tower
                '''
    #only need to return solution and len(solution)
    #return len(solution), solution
    parse_input = []
    for p in solution:
        parse_input.append(str(p.x + " " + p.y))
    solution = Solution.parse(parse_input, instance)
    return solution

    #From solution parse the coordinates of each tower

SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive
}


# You shouldn't need to modify anything below this line.
def infile(args):
    if args.input == "-":
        return StdinFileWrapper()

    return Path(args.input).open("r")


def outfile(args):
    if args.output == "-":
        return StdoutFileWrapper()

    return Path(args.output).open("w")


def main(args):
    with infile(args) as f:
        instance = Instance.parse(f.readlines())
        solver = SOLVERS[args.solver]
        solution = solver(instance)
        assert solution.valid()
        with outfile(args) as g:
            print("# Penalty: ", solution.penalty(), file=g)
            solution.serialize(g)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve a problem instance.")
    parser.add_argument("input", type=str, help="The input instance file to "
                        "read an instance from. Use - for stdin.")
    parser.add_argument("--solver", required=True, type=str,
                        help="The solver type.", choices=SOLVERS.keys())
    parser.add_argument("output", type=str,
                        help="The output file. Use - for stdout.",
                        default="-")
    main(parser.parse_args())
