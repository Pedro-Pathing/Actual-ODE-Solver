package MathUtil;

import org.apache.commons.math3.util.FastMath;

import java.awt.geom.Point2D;

public class SquareCoordinates {
    private double distance = 0;
    private double theta = 0;

    public SquareCoordinates(double distance, double theta) {
        this.distance = distance;
        this.theta = theta;
    }

    public SquareCoordinates(Point2D.Double cartesianPoint, boolean isRed) {
        distance = FastMath.min(FastMath.abs(72 - cartesianPoint.x), FastMath.abs(72 - cartesianPoint.y));
        if (!isRed) {
            theta = normalizeAngleBlue(FastMath.atan2(cartesianPoint.y,cartesianPoint.x));
        } else {
            theta = normalizeAngleBlue(FastMath.atan2(cartesianPoint.y,cartesianPoint.x));
        }
    }

    private double normalizeAngleBlue(double angleRadians) {
        double angle = angleRadians;
        while (angle < 0) angle += 2 * Math.PI;
        while (angle > Math.PI) angle -= 2 * Math.PI;
        return angle;
    }

    private double normalizeAngleRed(double angleRadians) {
        double angle = angleRadians;
        while (angle < Math.PI) angle += 2 * Math.PI;
        while (angle > 2 * Math.PI) angle -= 2 * Math.PI;
        return angle;
    }

    public double getX() {
        if (FastMath.cos(theta) >= FastMath.sin(theta)) {
            return 72 + distance/Math.tan(theta);
        } else {
            return 72 + distance;
        }
    }

    public double getY() {
        if (FastMath.sin(theta) >= FastMath.cos(theta)) {
            return 72 + distance * Math.tan(theta);
        } else {
            return 72 + distance;
        }
    }

    public double getDistance() {
        return distance;
    }

    public double getTheta() {
        return theta;
    }

    public void setDistance(double distance) {
        this.distance = distance;
    }

    public void setTheta(double theta) {
        this.theta = theta;
    }

    public Point2D.Double getCartesian() {
        return new Point2D.Double(getX(), getY());
    }
}
