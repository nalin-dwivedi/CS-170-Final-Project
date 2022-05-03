# Spring 2022 CS170 Project (Team 2: pay2win)

## Spec
https://cs170.org/assets/pdf/project_spec.pdf

## Requirements

The project was developed using Python 3.9, but it should work with Python versions 3.6+.

## How to View Outputs

Clone the repository onto your local computer and begin by deleting the `outputs/` subdirectory. **Note: you will have to delete the** `outputs/` **subdirectory every time you want to rerun our algorithm against the inputs in the** `inputs/` **subdirectory.** Below you will find procedures to produce our highest performing output files for small, medium and large inputs respectively.

### For Small and Medium Inputs

Navigate to the `python/solve_all.py` file and first uncomment the if-elif block from lines 38-41 in `def solver(size: Size, instance: Instance) -> Solution` that correspond to the if-elif statements `if size == Size.SMALL` and `elif size == Size.MEDIUM`. Make sure that all other lines in the `solver` function are commented out. Then change the return statements for both if statements on lines 39 and 41 to `return algo_ver2(instance)`. After making sure the `outputs/` subdirectory is deleted, open a terminal window and run `python3 python/solve_all.py inputs outputs` in the repo directory which will create `outputs/small/` and `outputs/medium/` folders filled with `.out` files corresponding to the `.in` files in the `inputs/small` and `inputs/medium` subdirectories. 

Please store these output files somewhere safe to view them later as on every rerun, the `outputs/` subdirectory will be overriden. Note the `outputs/large/` folder will be empty because we have not ran the large inputs using the `algo_ver2` algorithm. These outputs will be produced with the procedure below.

### For Large Inputs

Navigate to the `python/solve_all.py` file and first uncomment the if block from lines 42-43 in `def solver(size: Size, instance: Instance) -> Solution` that correspond to the if statement `if size == Size.LARGE`. Make sure that all other lines in the `solver` function are commented out. Then change the return statement for the if statement on line 43 to `return algo_ver3(instance)`. After making sure the `outputs/` subdirectory is deleted, open a terminal window and run `python3 python/solve_all.py inputs outputs` in the repo directory which will create an `outputs/large/` folder filled with `.out` files corresponding to the `.in` files in the `inputs/large/` subdirectory. 

Note the `outputs/small/` and `outputs/medium/` folders within the `outputs/` folder will be empty because we have not ran the small and medium inputs using the `algo_ver3` algorithm. These outputs will be produced with the procedure above.
