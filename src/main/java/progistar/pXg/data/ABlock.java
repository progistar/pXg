package progistar.pXg.data;

/**
 * Annotation Block <br>
 * 
 * @author gistar
 *
 */
public class ABlock implements Comparable<ABlock>{

	public int transcriptIndex;
	public int start;
	public int end;
	public byte feature;
	public boolean strand;
	
	@Override
	public int compareTo(ABlock o) {
		if(this.start < o.start) {
			return -1;
		}else if(this.start > o.start) {
			return 1;
		}
		
		return 0;
	}
}
