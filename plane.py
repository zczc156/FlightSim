import numpy as np


class Plane:
    def __init__(self, prop, aero_center, m, Iz, AR, A, T):
        """
        Parameters
        ----------
        prop : (float, float)
            Position of propeller (thurst) relative to center of gravity.
        aero_center : (float, float)
            Position of aerodynamic center relative to the center of gravity.
        m : float
            Mass of the plane.
        Iz : float
            Moment of inertia of the plane.
        AR : float
            Aspect ratio of the airfoil.
        A : float
            Area of the whold plane.
        T : float
            Thrust that is put to the plane.
        """
        self.xp, self.yp = prop
        self.xa, self.ya = aero_center
        self.m = m
        self.Iz = Iz
        self.AR = AR
        self.A = A
        self.T = T

    def lift_coef(self, alpha):
        """The function alculates lift coefficient of the plane using AOA alpha.

        Parameters
        ----------
        alpha : float
            Angle of attack (AOA) in rad.

        Returns
        -------
        float
            Lift coefficient of the plane.
        """
        return 2 * np.pi * np.sin(alpha) / (1 + 2 / self.AR)

    def drag_coef(self, alpha):
        """The function alculates drag coefficient of the plane using AOA alpha.

        Parameters
        ----------
        alpha : float
            Angle of attack (AOA) in degrees.

        Returns
        -------
        float
            Drag coefficient of the plane.
        """
        drag_coef_inf = 0.01
        lift_coef = self.lift_coef(alpha)

        return drag_coef_inf + lift_coef ** 2 / (np.pi * self.AR)


if __name__ == "__main__":
    plane = Plane((0.01, 0.01), (-0.01, 0), 0.1, 0.001, 3, 0.3, 4)

    cl = plane.lift_coef(np.deg2rad(5))
    cd = plane.drag_coef(np.deg2rad(5))
    print(cl, cd)
