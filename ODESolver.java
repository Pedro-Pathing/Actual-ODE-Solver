import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class ODESolver {

    public static class VelocityAndThetaODE implements FirstOrderDifferentialEquations {
        private double vmax;    // Maximum velocity in inches/second
        private double mu_k;    // Coefficient of kinetic friction
        private double N;       // Normal force in Newtons
        private double c1;      // Sinusoidal coefficient

        public VelocityAndThetaODE(double vmax, double mu_k, double N, double c1) {
            this.vmax = vmax;
            this.mu_k = mu_k;
            this.N = N;
            this.c1 = c1;
        }

        @Override
        public int getDimension() {
            return 2; // x and y components only
        }

        @Override
        public void computeDerivatives(double t, double[] state, double[] dState) {
            double x = state[0];
            double y = state[1];

            double frictionForce = mu_k * N + c1 * Math.abs(y);
            dState[0] = vmax - frictionForce; // dx/dt
            dState[1] = 0; // dy/dt (constant heading, no vertical movement)
        }
    }

    public static double calculateArcLength(List<double[]> waypoints) {
        double arcLength = 0.0;

        for (int i = 1; i < waypoints.size(); i++) {
            double[] p1 = waypoints.get(i - 1);
            double[] p2 = waypoints.get(i);

            double distance = Math.sqrt(Math.pow(p2[0] - p1[0], 2) + Math.pow(p2[1] - p1[1], 2));
            arcLength += distance;
        }

        return arcLength;
    }

    public static double[] calculateWaypointPosition(double P0, double P1, double P2, double P3, double s) {
        double x = Math.pow(s, 3) * (P3 - 3 * P2 + 3 * P1 - P0) +
                Math.pow(s, 2) * (3 * P2 - 6 * P1 + 3 * P0) +
                s * (3 * P1 - 3 * P0) +
                P0;
        double y = Math.pow(s, 2) * (P3 - 3 * P2 + 3 * P1 - P0) +
                6 * s * (P2 - 2 * P1 + P0) +
                (3 * P1 - P0);
        return new double[]{x, y};
    }

    public static boolean validateWaypoint(double x, double y, double a, double h, double W_s, double b, double k, double L_s) {
        return a <= x && x <= 144 - a &&
                a <= y && y <= 72 - a &&
                Math.abs(x - h) >= W_s / 2 + b / 2 &&
                Math.abs(y - k) >= L_s / 2 + b / 2;
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Inputs
        System.out.println("Enter maximum velocity (vmax) in inches/second: ");
        double vmax = scanner.nextDouble();

        System.out.println("Enter kinetic friction coefficient (mu_k): ");
        double mu_k = scanner.nextDouble();

        System.out.println("Enter normal force (N) in Newtons: ");
        double N = scanner.nextDouble();

        System.out.println("Enter sinusoidal coefficient (c1): ");
        double c1 = scanner.nextDouble();

        System.out.println("Enter initial heading angle (theta0) in degrees: ");
        double theta0Degrees = scanner.nextDouble();
        double theta0 = Math.toRadians(theta0Degrees);

        System.out.println("Enter target heading angle (H) in degrees: ");
        double HDegrees = scanner.nextDouble();
        double H = Math.toRadians(HDegrees);

        System.out.println("Enter angular velocity (ω) in radians/second: ");
        double omega = scanner.nextDouble();

        System.out.println("Enter waypoint parameters (a, h, W_s, b, k, L_s): ");
        double a = scanner.nextDouble();
        double h = scanner.nextDouble();
        double W_s = scanner.nextDouble();
        double b = scanner.nextDouble();
        double k = scanner.nextDouble();
        double L_s = scanner.nextDouble();

        System.out.println("Enter waypoint positions (P0, P1, P2, P3): ");
        double P0 = scanner.nextDouble();
        double P1 = scanner.nextDouble();
        double P2 = scanner.nextDouble();
        double P3 = scanner.nextDouble();

        // Step 1: Calculate dTheta and t2
        double dTheta = Math.abs(H - theta0);
        double t2 = dTheta / omega;
        System.out.println("Time to align heading (t2): " + t2 + " seconds");

        // Step 2: Generate and validate waypoints
        List<double[]> waypoints = new ArrayList<>();
        for (double s = 0; s <= 1; s += 0.01) { // s ranges from 0 to 1
            double[] waypoint = calculateWaypointPosition(P0, P1, P2, P3, s);
            if (validateWaypoint(waypoint[0], waypoint[1], a, h, W_s, b, k, L_s)) {
                waypoints.add(waypoint);
            } else {
                System.out.println("Waypoint validation failed at s = " + s);
            }
        }

        // Step 3: Calculate arc length (L)
        double L = calculateArcLength(waypoints);
        System.out.println("Arc length (L): " + L + " inches");

        // Step 4: Solve ODE for t1
        FirstOrderDifferentialEquations ode = new VelocityAndThetaODE(vmax, mu_k, N, c1);
        double[] initialState = {0.0, 0.0}; // Starting position
        double t0 = 0.0;
        double t1Guess = 10.0; // Initial guess for t1 (time to traverse path)

        FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-6, 1.0, 1.0e-10, 1.0e-10);
        integrator.integrate(ode, t0, initialState, t1Guess, initialState);

        // Step 5: Output results
        double t1 = t1Guess;
        System.out.println("Time to traverse path (t1): " + t1 + " seconds");
        System.out.println("Change in heading (dTheta): " + Math.toDegrees(dTheta) + " degrees");
        System.out.println("Rate: -dt1/dTheta = " + (-t1 / dTheta) + " seconds/radian");
    }
}