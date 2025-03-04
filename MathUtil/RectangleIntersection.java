package MathUtil;

import java.awt.geom.Point2D;

/**
 * Utility class for rectangle operations, including rotated rectangle creation and overlap checks.
 */
public class RectangleIntersection {

    public static class Rectangle {
        Point2D.Double[] vertices;

        public Rectangle(Point2D.Double[] vertices) {
            this.vertices = vertices;
        }
    }

    /**
     * Rotates a point around a pivot by a given angle (in degrees).
     *
     * @param p Point to rotate.
     * @param pivot Pivot point to rotate around.
     * @param angle Angle in degrees.
     * @return New rotated point.
     */
    public static Point2D.Double rotatePoint(Point2D.Double p, Point2D.Double pivot, double angle) {
        double rad = Math.toRadians(angle);
        double cos = Math.cos(rad);
        double sin = Math.sin(rad);

        double translatedX = p.x - pivot.x;
        double translatedY = p.y - pivot.y;

        double rotatedX = translatedX * cos - translatedY * sin;
        double rotatedY = translatedX * sin + translatedY * cos;

        return new Point2D.Double(rotatedX + pivot.x, rotatedY + pivot.y);
    }

    /**
     * Creates a rectangle centered at a given point with specified width, height, and rotation angle.
     *
     * @param center Center point of the rectangle.
     * @param width Width of the rectangle.
     * @param height Height of the rectangle.
     * @param rotation Rotation angle in degrees.
     * @return A rectangle object with vertices computed based on the rotation.
     */
    public static Rectangle createRotatedRectangle(Point2D.Double center, double width, double height, double rotation) {
        double halfWidth = width / 2;
        double halfHeight = height / 2;

        Point2D.Double[] vertices = new Point2D.Double[4];
        vertices[0] = new Point2D.Double(center.x - halfWidth, center.y - halfHeight);
        vertices[1] = new Point2D.Double(center.x + halfWidth, center.y - halfHeight);
        vertices[2] = new Point2D.Double(center.x + halfWidth, center.y + halfHeight);
        vertices[3] = new Point2D.Double(center.x - halfWidth, center.y + halfHeight);

        for (int i = 0; i < vertices.length; i++) {
            vertices[i] = rotatePoint(vertices[i], center, rotation);
        }

        return new Rectangle(vertices);
    }

    /**
     * Computes the maximum overlap between two rectangles using the Separating Axis Theorem (SAT).
     * If the rectangles do not overlap, returns 0.
     *
     * @param rect1 First rectangle.
     * @param rect2 Second rectangle.
     * @return Maximum overlap distance if rectangles overlap, or 0 if they are separated.
     */
    public static double getMaximumOverlap(Rectangle rect1, Rectangle rect2) {
        double maxOverlap = 0; // Start with 0, will maximize across all axes.

        // Check all edges of rect1
        for (int i = 0; i < rect1.vertices.length; i++) {
            Point2D.Double p1 = rect1.vertices[i];
            Point2D.Double p2 = rect1.vertices[(i + 1) % rect1.vertices.length];

            Point2D.Double axis = computePerpendicularAxis(p1, p2);

            if (!checkAxisOverlap(axis, rect1, rect2)) {
                return 0;  // Separating axis found - no overlap
            }

            double overlap = computeAxisOverlap(axis, rect1, rect2);
            maxOverlap = Math.max(maxOverlap, overlap);
        }

        // Check all edges of rect2
        for (int i = 0; i < rect2.vertices.length; i++) {
            Point2D.Double p1 = rect2.vertices[i];
            Point2D.Double p2 = rect2.vertices[(i + 1) % rect2.vertices.length];

            Point2D.Double axis = computePerpendicularAxis(p1, p2);

            if (!checkAxisOverlap(axis, rect1, rect2)) {
                return 0;  // Separating axis found - no overlap
            }

            double overlap = computeAxisOverlap(axis, rect1, rect2);
            maxOverlap = Math.max(maxOverlap, overlap);
        }

        return maxOverlap;  // Maximum overlap across all axes.
    }

    private static Point2D.Double computePerpendicularAxis(Point2D.Double p1, Point2D.Double p2) {
        Point2D.Double edge = new Point2D.Double(p2.x - p1.x, p2.y - p1.y);
        return new Point2D.Double(-edge.y, edge.x);
    }

    private static boolean checkAxisOverlap(Point2D.Double axis, Rectangle rect1, Rectangle rect2) {
        double[] proj1 = projectRectangle(rect1, axis);
        double[] proj2 = projectRectangle(rect2, axis);
        return proj1[1] >= proj2[0] && proj2[1] >= proj1[0];  // True if projections overlap
    }

    private static double computeAxisOverlap(Point2D.Double axis, Rectangle rect1, Rectangle rect2) {
        double[] proj1 = projectRectangle(rect1, axis);
        double[] proj2 = projectRectangle(rect2, axis);
        return Math.max(0, Math.min(proj1[1], proj2[1]) - Math.max(proj1[0], proj2[0]));
    }

    /**
     * Projects a rectangle onto a given axis and returns the min/max projection values.
     *
     * @param rect Rectangle to project.
     * @param axis Axis vector.
     * @return Array containing {minProjection, maxProjection}.
     */
    private static double[] projectRectangle(Rectangle rect, Point2D.Double axis) {
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;

        double axisLength = Math.sqrt(axis.x * axis.x + axis.y * axis.y);

        for (Point2D.Double vertex : rect.vertices) {
            double projection = (vertex.x * axis.x + vertex.y * axis.y) / axisLength;
            min = Math.min(min, projection);
            max = Math.max(max, projection);
        }

        return new double[]{min, max};
    }
}