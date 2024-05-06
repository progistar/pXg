package progistar.pXg.data;

import progistar.pXg.constants.Constants;
import progistar.pXg.utils.IndexConvertor;

public class Mutation {

	public byte type;

	public String altSeq;
	public String refSeq;
	public int relPos;

	public int genomicPosition;
	public int chrIndex;

	@Override
	public String toString () {
		if(type == Constants.SNP) {
			return IndexConvertor.indexToChr(chrIndex) +":" +genomicPosition+refSeq+">"+altSeq;
		} else if(type == Constants.INS){
			return IndexConvertor.indexToChr(chrIndex) +":" +genomicPosition+"ins"+altSeq;
		} else if(type == Constants.DEL){
			return IndexConvertor.indexToChr(chrIndex) +":" +genomicPosition+"del"+refSeq;
		}

		return "NA";
	}
}
