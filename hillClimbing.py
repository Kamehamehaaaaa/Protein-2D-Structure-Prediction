import random
import matplotlib.pyplot as plt
import time

hpSequence = "HPHPHPPPPHPPHPPH"

moves = ['L', 'F', 'R']

def generateGrid(path):
    coords = [[0,0]]
    prevTile = (-1,0) # left-to-right start

    for i in range(len(path)):
        new_coords = []
        currTile = coords[-1]
        try:
            match path[i]:
                case 'L':
                    # whats not changing will change
                    if prevTile[0] - currTile[0] == 0:
                        x = prevTile[1] - currTile[1] + currTile[0]
                        y = currTile[1]
                    elif prevTile[1] - currTile[1] == 0:
                        x = currTile[0]
                        y = currTile[0] - prevTile[0] + currTile[1]
                    else:
                        raise ValueError("Error with left move")
                case 'R':
                    # whats not changing will change
                    if prevTile[0] - currTile[0] == 0:
                        x = currTile[1] - prevTile[1] + currTile[0]
                        y = currTile[1]
                    elif prevTile[1] - currTile[1] == 0:
                        x = currTile[0]
                        y = prevTile[0] - currTile[0] + currTile[1]
                    else:
                        raise ValueError("Error with right move")
                case 'F':
                    # whats same will be same 
                    x = currTile[0] - prevTile[0] + currTile[0]
                    y = currTile[1] - prevTile[1] + currTile[1]
                case _:
                    raise ValueError("Invalid Move!!")
        except ValueError as e:
            print("Invalid move: ", e)
            exit()

        new_coords.append(x)
        new_coords.append(y)

        # check the Self-Avoiding Walk constraint
        if new_coords in coords:
            return False, []
        coords.append(new_coords)
        prevTile = currTile

    return True, coords

def calculateHHContact(sequence, grid):
    energy = 0
    directions = [(1,0), (0,1), (-1,0), (0,-1)]

    for i in range(len(sequence)):
        if sequence[i] == 'H':
            pos = grid[i]
            neighbours = [[pos[0]+j[0], pos[1]+j[1]] for j in directions]

            for j in neighbours:
                if j in grid:
                    idx = grid.index(j)
                    if sequence[idx] == 'H':
                        energy += 1

            if i != 0 and sequence[i-1] == 'H':
                energy -= 1
            if i != len(sequence) - 1 and sequence[i+1] == 'H':
                energy -= 1
    
    # since energy is going to be counted twice
    return -1 * (energy//2)

def randomPath(length):
    while True:
        path = [random.choice(moves) for _ in range(length)]
        # path = []
        # for i in "LRFRFRFLRFFFRRL":
        #     path.append(i)
        valid, grid = generateGrid(path)
        if valid:
            return path, grid
            

# _, grid = generateGrid("LRRFFLLF")
# seq = "HPHHHPHHP"
# energy = calculateHHContact(seq, grid)

# print(grid)
# print(energy)

# path, grid = randomPath(len(seq)-1)
# energy = calculateHHContact(seq, grid)

# print(grid)
# print(path, energy)

def hill_climb(seq, max_iteration=1000):
    start = time.time()
    path, grid = randomPath(len(seq)-1)
    print("Random Fold: ", ''.join(path))
    bestEnergy = calculateHHContact(seq, grid)
    # print(path, bestEnergy)
    no_improvement = 0

    for iter in range(max_iteration):
        i = random.randrange(len(path))
        new_path = path[:]
        new_path[i] = random.choice([m for m in moves if m != path[i]])
        
        valid, new_grid = generateGrid(new_path)
        if not valid:
            continue
        new_energy = calculateHHContact(seq, new_grid)
        if new_energy < bestEnergy:
            bestEnergy = new_energy
            path = new_path[:]
            grid = new_grid[:]
            no_improvement = 0
        else:
            no_improvement += 1
        # premature termination when no improvement for 100 continuous cycles
        if no_improvement == 100:
            break
        
        # print(path, new_energy, bestEnergy)
    
    
    return bestEnergy, path, grid, iter, time.time()-start


def getHPLattice(sequence, grid):
    plt.figure(figsize=(6,6))

    # Plot bonds (edges)
    for i in range(len(grid)-1):
        x = [grid[i][0], grid[i+1][0]]
        y = [grid[i][1], grid[i+1][1]]
        plt.plot(x, y, color='gray', linewidth=2, zorder=1)

    # Plot residues (nodes)
    for i, (x, y) in enumerate(grid):
        if sequence[i] == 'H':
            plt.scatter(x, y, color='orange', s=200, edgecolors='black', zorder=2, label='Hydrophobic' if i==0 else "")
        else:
            plt.scatter(x, y, color='skyblue', s=200, edgecolors='black', zorder=2, label='Polar' if i==0 else "")
        plt.text(x, y+0.2, str(i), fontsize=8, ha='center')

    # Make lattice grid look clean
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.gca().set_aspect('equal')
    plt.title("2D Lattice Protein Fold Visualization")
    plt.legend(loc='upper left')
    plt.xlabel("x-axis")
    plt.ylabel("y-axis")
    plt.tight_layout()
    plt.show()


# print(hill_climb(seq, 10))

# seq = "HPHPPHHPHPPHPHHPPHPH"
# seq = "HPHPPHHPHPPH"
# seq = "HHPPHPHPH"
# seq = "HPPHHPHPHHPPHPHPPHHPPHHPHPHPHHPPPHHPPHPHPHHPPHPPHPHHPPHPPHHPPHPHPPHPPHPHHPPPHHPPHPPHPPHPHHPPHPH"
seq = "HHHPPPPHPHPHPPHH"

energies = []
iters = []
best_energy = 0
best_grids = []
iters = []
runtimes = []
for simulation in range(30):
    e,path,grid, iter, runtime= hill_climb(seq)
    if e < best_energy:
        best_energy = e
        best_grids = [grid]
    elif e == best_energy:
        best_grids.append(grid)
    energies.append(e)
    iters.append(iter)
    runtimes.append(runtime)
    print(f"Iteration {simulation}:\nFinal Fold : {''.join(path)}\nFinal Energy: {e}")

print(f"Average energy: {sum(energies)/30}")
print(f"Average iters: {sum(iters)/30}")
print(f"Average Runtime: {sum(runtimes)/30}")

plt.plot(energies, marker='o')
plt.title("Hill Climbing: Protein Folding Energies")
plt.xlabel("Trial")
plt.ylabel("Energy")
plt.show()

print(best_grids)
print(iters)
for grid in best_grids:
    getHPLattice(seq, grid)

            
            
        
    
