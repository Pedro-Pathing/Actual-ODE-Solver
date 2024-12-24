import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class ODESolver {

    // Compute the velocity vector based on the provided velocity equation
    public static double[] computeVelocity(double[] PPrime, double[] R, double theta, double vmax, double mu_k, double N, double c1) {
        double cosTheta = FastMath.cos(theta);
        double sinTheta = FastMath.sin(theta);
        double magnitude = vmax - (PPrime[0] / R[0]) * mu_k * N + c1 * Math.abs(sinTheta);

        return new double[]{PPrime[0] * magnitude * cosTheta, PPrime[1] * magnitude * sinTheta};
    }

    // Compute waypoint using Lagrange interpolation
    public static double[] computeWaypoint(double[] x, double[] y, double t) {
        double xt = 0, yt = 0;
        for (int i = 0; i < x.length; i++) {
            double li = 1;
            for (int j = 0; j < x.length; j++) {
                if (i != j) {
                    li *= (t - j) / (i - j);
                }
            }
            xt += li * x[i];
            yt += li * y[i];
        }
        return new double[]{xt, yt};
    }

    // Validate a waypoint using the provided constraints
    public static boolean validateWaypoint(double x, double y, double a, double h, double W_s, double b, double k, double L_s) {
        return a <= y && y <= 144 - a &&
                a <= x && x <= 72 - a &&
                Math.abs(y - h) >= W_s / 2 + b / 2 &&
                Math.abs(x - k) >= L_s / 2 + b / 2;
    }

    // Function to compute dt2/dTheta
    public static double computeDt2OverDTheta(double theta0, double theta1, double omega) {
        double dTheta = Math.abs(theta1 - theta0);
        return dTheta / omega;
    }

    // Function for -dt1/dTheta
    public static double calculateDt1OverDTheta(double theta, double vmax, double mu_k, double N, double c1) {
        return -1 / (vmax - mu_k * N + c1 * Math.abs(Math.sin(theta)));
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // User inputs
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

        System.out.println("Enter waypoint x-coordinates (comma or space-separated): ");
        scanner.nextLine(); // Clear buffer
        String[] xCoordsInput = scanner.nextLine().split("[,\\s]+");
        double[] xCoords = new double[xCoordsInput.length];
        for (int i = 0; i < xCoordsInput.length; i++) {
            xCoords[i] = Double.parseDouble(xCoordsInput[i]);
        }

        System.out.println("Enter waypoint y-coordinates (comma or space-separated): ");
        String[] yCoordsInput = scanner.nextLine().split("[,\\s]+");
        double[] yCoords = new double[yCoordsInput.length];
        for (int i = 0; i < yCoordsInput.length; i++) {
            yCoords[i] = Double.parseDouble(yCoordsInput[i]);
        }

        // Calculate dt2/dTheta
        double dt2OverDTheta = computeDt2OverDTheta(theta0, H, omega);

        // Generate and validate waypoints
        List<double[]> validWaypoints = new ArrayList<>();
        for (double t = 0; t <= 1; t += 0.01) { // t ranges from 0 to 1
            double[] waypoint = computeWaypoint(xCoords, yCoords, t);
            if (validateWaypoint(waypoint[0], waypoint[1], a, h, W_s, b, k, L_s)) {
                validWaypoints.add(waypoint);
            }
        }

        if (validWaypoints.isEmpty()) {
            System.out.println("No valid waypoints found!");
            return;
        }

        double[] P1 = validWaypoints.get(1);
        double[] P2 = validWaypoints.get(validWaypoints.size() - 2);

        System.out.println("First valid waypoint (P1): (" + P1[0] + ", " + P1[1] + ")");
        System.out.println("Last valid waypoint (P2): (" + P2[0] + ", " + P2[1] + ")");

        // Calculate arc length (L) using Simpson's integrator
        SimpsonIntegrator integrator = new SimpsonIntegrator();
        double L = integrator.integrate(10000, x -> {
            double[] velocity = computeVelocity(new double[]{1, 0}, new double[]{1, 0}, theta0, vmax, mu_k, N, c1);
            return Math.sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1]);
        }, 0, 1);

        System.out.println("Arc length (L): " + L + " inches");

        // Calculate total time t1
        double t1 = L / vmax;
        System.out.println("Time to traverse path (t1): " + t1 + " seconds");

        // Output the rate function
        System.out.println("Rate (dt1/dTheta = -dt2/dTheta): ");
        for (double theta = theta0; theta <= H; theta += 0.1) {
            double rate = -dt2OverDTheta;
            System.out.println("Theta (radians): " + theta + ", Rate: " + rate + " seconds/radian");
        }

        // Total time
        double totalTime = t1 + (dt2OverDTheta * Math.abs(H - theta0));
        System.out.println("Total time: " + totalTime + " seconds");
    }
}