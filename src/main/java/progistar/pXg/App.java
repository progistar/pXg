package progistar.pXg;

import progistar.pXg.constants.Parameters;
import progistar.pXg.data.parser.ParameterParser;
import progistar.pXg.processor.Master;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
    	ParameterParser.parseParams("params.txt");
    	
    	Master.ready(Parameters.genomicAnnotationFilePath, Parameters.sequenceFilePath, Parameters.peptideFilePath);
    	Master.run();
    }
}
