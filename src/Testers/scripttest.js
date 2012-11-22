/* This file is read by scripttest.java
 *
 * Symbols 'spam' and 'eggs' were 'put' into engine prior to 'eval'
 */

function getNewName (s) {
    var out = s.replace("Spam", "Rotten Eggs");
    out = out.replace("My", "Your");
    return out;
}

function getNewAge (a) {
    return a / 3 + 4;
}

var new_name = getNewName(spam.get("name"));
var new_age = getNewAge(spam.get("age"));

eggs.put("name", new_name);
eggs.put("age", new_age);
