import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;


class PlotLineChart extends ApplicationFrame{
    DefaultCategoryDataset dataset = new DefaultCategoryDataset();
    JFreeChart lineChart;

    PlotLineChart( String applicationTitle , String chartTitle ) {
        super(applicationTitle);
    }

    void setFinalData(String chartTitle){
        JFreeChart lineChart = ChartFactory.createLineChart(
                chartTitle,
                "Rauschen","Genauigkeit",
                dataset,
                PlotOrientation.VERTICAL,
                true,true,false);

        ChartPanel chartPanel = new ChartPanel( lineChart );
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        setContentPane( chartPanel );
    }

    void addToDataset(double data, String name, String columnKey){
        dataset.addValue(data, name, columnKey);
    }

    JFreeChart getLineCHart(){
        return lineChart;
    }

    DefaultCategoryDataset getDataset(){
        return dataset;
    }
}
