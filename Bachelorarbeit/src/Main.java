import org.apache.commons.math3.stat.regression.GLSMultipleLinearRegression;

import java.util.ArrayList;
import java.util.List;

public class Main {

    /**gibt true zurueck, falls der Abstand des PositionPairs zu dem Punkt mit den Koordinaten positionX und positionY
     * kleiner als border ist.
     * Ansonsten true
     *
     * @param positionPairs
     * @param positionX
     * @param positionY
     * @param border
     * @return
     */
    private static boolean tooClose(ArrayList<NewPair> positionPairs, int positionX, int positionY, int border){
        for (NewPair p : positionPairs) {
            double pythagoras = Math.sqrt((Math.pow(positionX - p.getFirst(), 2) + (Math.pow(positionY - p.getSecond(), 2))));
            if(pythagoras < border){
                return true;
            }
        }
        return false;
    }

    /**Erzeugt eine Zufallszahl um den festgelegten Bereich.
     *
     * @param streuungsreichweite
     * @param bereich
     * @return
     */
    private static int getRandomPosition(int streuungsreichweite, int bereich){
        return (int)((streuungsreichweite/2) + Math.random() * bereich);
    }

    public static void main(String[] args) {
        List<List<NewPair>> clusters = new ArrayList<>();
        List<NewPair> c1 = new ArrayList<>();
        List<NewPair> c2 = new ArrayList<>();
        List<NewPair> c3 = new ArrayList<>();
        List<NewPair> c4 = new ArrayList<>();

        clusters.add(c1);
        clusters.add(c2);
        clusters.add(c3);
        clusters.add(c4);


        int streuungsreichweite = 15;
        int bereich = 75;
        int xK = getRandomPosition(streuungsreichweite, bereich);
        int yK = getRandomPosition(streuungsreichweite, bereich);
        ArrayList<NewPair> clusterPositions = new ArrayList<>();
        int first;
        int second;
        for (List<NewPair> x : clusters) {

            while(tooClose(clusterPositions, xK, yK, 2 * streuungsreichweite)){
                xK = getRandomPosition(streuungsreichweite, bereich);
                yK = getRandomPosition(streuungsreichweite, bereich);
            }

            for(int i = 0; i < 100; i++){
                first = (int)((xK - (streuungsreichweite/2)) + (Math.random() * streuungsreichweite));
                second = (int)((yK - (streuungsreichweite/2)) + (Math.random() * streuungsreichweite));
                NewPair p = new NewPair(first, second);
                x.add(p);
            }
            clusterPositions.add(new NewPair(xK,yK));
        }

        /*int k = 0;
        for (List<NewPair> x : clusters) {
            System.out.println("Cluster" + ++k);
            for (NewPair p : x) {
                System.out.println(p.getFirst() + "," + p.getSecond());
            }
        }*/

        new DrawPoints(clusters);

        GLSMultipleLinearRegression regression = new GLSMultipleLinearRegression();
    }

}
