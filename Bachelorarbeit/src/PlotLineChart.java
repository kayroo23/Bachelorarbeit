import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;


class PlotLineChart extends ApplicationFrame{
    DefaultCategoryDataset dataset = new DefaultCategoryDataset();

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

    private DefaultCategoryDataset createDataset( ) {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset( );
        dataset.addValue( 15 , "schools" , "1970" );
        dataset.addValue( 30 , "schools" , "1980" );
        dataset.addValue( 60 , "schools" ,  "1990" );
        dataset.addValue( 120 , "schools" , "2000" );
        dataset.addValue( 240 , "schools" , "2010" );
        dataset.addValue( 300 , "schools" , "2014" );
        return dataset;
    }

    void addToDataset(double data, String name, String columnKey){
        dataset.addValue(data, name, columnKey);
    }
}
