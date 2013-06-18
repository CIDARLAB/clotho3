/**
 * This JavaScript file is the initiation script for the Mind's ScriptEngine (the engine that executes submit commands).
 * 
 * It deals with String conversions during communication with the ssAPI and adds some sloppiness and responsiveness.
 * 
 * Only commands that enter the ssAPI via 'submit' ever go through these functions.
 ***/

var clotho = {};

function resolveToString(input) {
    if(typeof input == "string" ) {
        return input;
    } else {
        return JSON.stringify(input);
    }
}


//var complete = clotho.autocomplete("walk the");
//clotho.say(JSON.stringify(complete));
clotho.autocomplete = function(inputPhrase, callback) {
    var args = resolveToString(inputPhrase);
    var sresult = clothoJava.autocomplete(args);
    var result = JSON.parse(sresult);
    if(callback) callback(result );
    else return result;
};

clotho.autocompleteDetail = function(inputPhrase, callback) {
    var args = resolveToString(inputPhrase);
    var sresult = clothoJava.autocompleteDetail(args);
    var result = JSON.parse(sresult);
    if(callback) callback(result );
    else return result;
};

clotho.submit = function(inputPhrase, callback) {
    var args = resolveToString(inputPhrase);
    clothoJava.submit(args);
    if(callback) callback( );
};

clotho.clear = function(callback) {
    clothoJava.clear();
    if(callback) callback( );
}

clotho.learn = function(nativeCmd, formalCmd, callback) {
    var args0 = resolveToString(nativeCmd);
    var args1 = resolveToString(formalCmd);
    clothoJava.learn(nativeCmd, formalCmd);
    if(callback) callback( );
};

//clotho.say("Says come up in the activity log");
clotho.say = function(inputPhrase, callback) {
    var args = resolveToString(inputPhrase);
    clothoJava.say(args);
    if(callback) callback();
};

//function say(args, callback) {
//    clothoJava.say(args);
//    if(callback) callback();
//}

clotho.alert = function(inputPhrase, callback) {
    var args = resolveToString(inputPhrase);
    clothoJava.alert(args);
    if(callback) callback();
};

clotho.log = function(inputPhrase, callback) {
    var args = resolveToString(inputPhrase);
    clothoJava.log(args);
    if(callback) callback();
};

clotho.note = function(inputPhrase, callback) {
    var args = resolveToString(inputPhrase);
    clothoJava.note(args);
    if(callback) callback();
};

//var paris = clotho.get('51b8b90050765ba45b5a8217');
//clotho.say(paris.state);
clotho.get = function(sharableRef, callback) {
    var args = resolveToString(sharableRef);
    var sresult = clothoJava.get(args);

    var result = JSON.parse(sresult);

    if(callback) callback(result);
    else return result;
};

clotho.set = function(sharableData, callback) {
    var args = resolveToString(sharableData);
    var sresult = clothoJava.set(args);

    var result = JSON.parse(sresult);

    if(callback) callback(result);
    else return result;
};

//var inst = {"city":"Paris","className":"org.clothocad.model.Institution","country":"United States of America","isDeleted":false,"lastModified":{"$date":"2013-06-08T02:48:10.254Z"},"name":"Test institution","state":"Massachusetts"};clotho.say(inst.city);   var result = clotho.create(inst);
clotho.create = function(sharableData, callback) {
    var args = resolveToString(sharableData);
    var sresult = clothoJava.create(args);
    var result = JSON.parse(sresult);
    if(callback) callback(result );
    else return result;
};

clotho.destroy = function(sharableRef, callback) {
    var args = resolveToString(sharableRef);
    clothoJava.destroy(args);
    if(callback) callback( );
};

//clotho.query('city', 'Townsville');
clotho.query = function(query, callback) {
    var jsonstr = resolveToString(query);
    var slisty = clothoJava.query(jsonstr);
    //println('query has slisty as ' + slisty);
    var listy = JSON.parse(slisty);
    if(callback) callback(listy);
    else return listy;
};

clotho.edit = function(sharableRef, callback) {
    var args = resolveToString(sharableRef);
    clothoJava.edit(args);
    if(callback) callback();
};

clotho.run = function(funcRef, sharableData, callback) {
    var args0 = resolveToString(funcRef);
    var args1 = resolveToString(sharableData);
    clothoJava.run(args0, args1);
    if(callback) callback();
};

clotho.show = function(viewRef, sharableData, position, recipients) {
    var args0 = resolveToString(viewRef);
    var args1 = resolveToString(sharableData);
    var args2 = resolveToString(position);
    var args3 = resolveToString(recipients);
    clothoJava.show(args0, args1, args2, args3);
    if(callback) callback();
};

clotho.startTrail = function(trailRef) {
    clothoJava.startTrail(trailRef);
    if(callback) callback();
};
