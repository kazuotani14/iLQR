import csv
import numpy as np
import matplotlib.pyplot as plt 

def read_trajectory(n_states, n_controls, results_path):
    # TODO clean this up
    reader = csv.reader(open(results_path))        
    states = []
    controls = []
    for i, row in enumerate(reader):
        if i==0: continue    
        s = [float(num) for num in row[:n_states]]
        states.append(s)

        if row[-1] != ' ': # check for terminal state
            a = [float(num) for num in row[-n_controls:]]
            controls.append(a)

    states = np.array(states)
    controls = np.array(controls)
    return states, controls


if __name__ == "__main__":
    # TODO take n_states, n_controls, result_path from argparser

    # Currently set up for acrobot
    n_states = 4
    n_controls = 1

    results_path = 'build/ilqr_result.csv'
    x, u = read_trajectory(n_states, n_controls, results_path)
    plt.plot(x[:, 0], 'g', label='x1')
    plt.plot(x[:, 1], 'b', label='x2')
    plt.plot(u[:, 0], 'r--', label='u1')
    plt.legend()
    plt.show()

    # results_path = 'results/acrobot_sqr.csv'
    # x_sqr, a_sqr = read_trajectory(n_states, n_controls, results_path)
    # plt.subplot(121)
    # plt.plot(x_sqr[:, 0], 'g', label='x1')
    # plt.plot(x_sqr[:, 1], 'b', label='x2')
    # plt.plot(a_sqr[:, 0], 'r--', label='u1')
    # plt.title('Square cost function')
    # plt.legend()

    # results_path = 'results/acrobot_sabs.csv'
    # x_sabs, a_sabs = read_trajectory(n_states, n_controls, results_path)
    # plt.subplot(122)
    # plt.plot(x_sabs[:, 0], 'g', label='x1')
    # plt.plot(x_sabs[:, 1], 'b', label='x2')
    # plt.plot(a_sabs[:, 0], 'r--', label='u1')
    # plt.title('Smooth abs cost function')
    # plt.legend()

    # plt.savefig('acrobot_cost_comparison.png')
    # plt.show()

