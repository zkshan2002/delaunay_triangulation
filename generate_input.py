import numpy as np

if __name__ == '__main__':
    all_points = np.random.uniform(-10.0, 10.0, size=(1000, 3))
    all_points[:, 2] = 0
    with open('data/input1000.xyz', 'w') as f:
        np.savetxt(f, all_points, fmt='%.10lf')
