

import javax.validation.constraints.Pattern;
import lombok.Getter;
import lombok.Setter;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
public class SequenceSchema extends ObjBase {
    
    public SequenceSchema(String name, String sequence){
        super(name);
        this.sequence = sequence;
    }
    @Getter
    @Setter
    @Pattern(regexp="[ATUCGRYKMSWBDHVN]*", flags={Pattern.Flag.CASE_INSENSITIVE})
    String sequence;
}

