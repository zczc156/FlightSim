import numpy as np
import matplotlib.pyplot as plt


from plane import Plane


class FlightSim:
    def __init__(self, plane, t, delta_t, air_density=1.275):
        """
        All angles are in rad.

        Parameters
        ----------
        plane : Plane
        t : float
            Determines how long simulation calculates.
        delta_t : float
            Time step of calculation.
        """
        self.plane = plane
        self.t = t
        self.delta_t = delta_t
        self.air_density = air_density

        # position, velocity
        # [rx, ry, theta, vx, vy, theta_dot]
        self.state = [0, 1, np.deg2rad(10), 5, 0, 0]
        self.states = None

        self.outdir = "./output"

        '''# acceleration
        # [ax, ay, theta_ddot]
        self.output = [0, 0, 0]'''

    def set_init_state(self, init_state):
        self.state = init_state

    def simulate(self):
        step_num = int(self.t / self.delta_t)

        self.states = np.zeros((step_num, 6))
        for step in range(step_num):
            self.states[step] = self.state

            rx, ry, theta, vx, vy, theta_dot = self.state
            phi = np.math.atan2(vy, vx)
            alpha = theta - phi

            cl = self.plane.lift_coef(alpha)
            cd = self.plane.drag_coef(alpha)
            v_sq = vx ** 2 + vy ** 2

            L = 0.5 * cl * self.air_density * np.cos(alpha) * self.plane.A * v_sq
            D = 0.5 * cd * self.air_density * np.sin(alpha) * self.plane.A * v_sq

            xp, yp, xa, ya = np.abs([self.plane.xp, self.plane.yp, self.plane.xa, self.plane.ya])

            ax = (self.plane.T * np.cos(theta) + L * np.sin(phi) - D * np.cos(phi)) / self.plane.m
            ay = (self.plane.T * np.sin(theta) + L * np.cos(phi) + D * np.sin(phi)) / self.plane.m - 9.81
            theta_ddot = (self.plane.T * np.sin(theta) * xp - self.plane.T * np.cos(theta) * yp
                          - L * np.cos(phi) * xa + L * np.sin(phi) * ya
                          + D * np.sin(phi) * xa - D * np.cos(phi) * ya) / self.plane.Iz

            vx += self.delta_t * ax
            vy += self.delta_t * ay
            theta_dot += self.delta_t * theta_ddot

            rx += self.delta_t * vx
            ry += self.delta_t * vy
            theta += self.delta_t * theta_dot

            self.state = [rx, ry, theta, vx, vy, theta_dot]

        return self.states

    def save_result(self, filename):
        with open(self.outdir + "/" + filename, "w") as f:
            lines_num = self.states.shape[0]
            for i in range(lines_num):
                line = "%.4f %.4f %.4f\n" % (self.states[i, 0], self.states[i, 1], self.states[i, 2])
                f.write(line)

    def plot_result(self, filename):
        states = []
        with open(self.outdir + "/" + filename, "r") as f:
            for i, line in enumerate(f):
                vals = line.split(' ')
                rx, ry, theta = float(vals[0]), float(vals[1]), float(vals[2])
                states.append([rx, ry, theta])
        states = np.array(states)

        maxval = np.amax(states, axis=0)
        minval = np.amin(states, axis=0)

        arr_len = 0.1
        for i in range(states.shape[0]):
            plt.arrow(states[i, 0], states[i, 1], arr_len * np.cos(states[i, 2]), arr_len * np.sin(states[i, 2]))
        plt.xlim(minval[0], maxval[0])
        plt.ylim(minval[1], maxval[1])
        plt.show()


if __name__ == "__main__":
    plane = Plane((0.1, -0.01), (-0.1, 0), 0.1, 0.001, 10, 0.5, 0)
    sim = FlightSim(plane, 1, 0.001)
    sim.simulate()
    sim.save_result("ex.txt")
    sim.plot_result("ex.txt")
