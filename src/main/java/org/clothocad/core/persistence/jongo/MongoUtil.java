package org.clothocad.core.persistence.jongo;

import com.google.common.collect.ImmutableMap;
import java.util.Map;

/**
 *
 * @author spaige
 */
public class MongoUtil {

    public static String sanitizeFieldName(String s) {
        String sanitized = s;
        for (String forbiddenCharacter : characterReplacements.keySet()) {
            sanitized = sanitized.replace(forbiddenCharacter, characterReplacements.get(forbiddenCharacter));
        }

        return sanitized;
    }

    public static String desanitizeFieldName(String s) {
        String desanitized = s;
        for (String forbiddenCharacter : characterReplacements.keySet()) {
            desanitized = desanitized.replace(characterReplacements.get(forbiddenCharacter), forbiddenCharacter);
        }

        return desanitized;
    }
    public static final Map<String, String> characterReplacements = ImmutableMap.of(".", "\uff0e", "$", "\uff04");
    
}
