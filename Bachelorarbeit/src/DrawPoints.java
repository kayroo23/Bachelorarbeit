
import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.util.List;

public class DrawPoints extends JFrame {
    final AffineTransform FLIP_X_COORDINATE;
    private final List<List<NewPair>> clusters;
    private final int scalingFactor = 5;


    public DrawPoints(List<List<NewPair>> c) {
        super("test");
        clusters = c;
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setSize(1000, 750);
        //setzt Koordinatenurspung nach unten links
        FLIP_X_COORDINATE = new AffineTransform(1, 0, 0, -1, 0, getHeight());
        setVisible(true);
    }


    @Override
    public void paint(Graphics graphic) {
        final int circleRadius = 3;

        Graphics2D g = (Graphics2D) graphic;
        g.setTransform(FLIP_X_COORDINATE);
        for (List<NewPair> x : clusters) {
            for (NewPair p : x) {
                g.drawOval(scalingFactor * p.getFirst(), scalingFactor * p.getSecond(), circleRadius,circleRadius);
            }
        }
    }


}