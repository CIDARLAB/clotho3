/**
 * This JavaScript file is the initiation script for the Mind's ScriptEngine (the engine that executes submit commands).
 * 
 * It deals with String conversions during communication with the ssAPI and adds some sloppiness and responsiveness.
 * 
 * Only commands that enter the ssAPI via 'submit' ever go through these functions.
 */
/*
Test sequence for testing CRUD operations:

//Create an object like Stephanie's test institution object
var newobj = clotho.create( {"lastModified":{"$date":"2013-06-13T02:47:51.288Z"},"name":"UCB","state":"C","className":"org.clothocad.model.Institution","isDeleted":false,"country":"United States of America","city":"Townsville"} );

//Get the newobj by a get call (with different sloppiness)
clotho.get(newobj.id);

clotho.get(newobj);

var wrapper = {};
wrapper.data = newobj;
clotho.get(wrapper);

wrapper = {};
wrapper.data = newobj.id;
clotho.get(wrapper);

wrapper = [];
wrapper[0] = newobj;
clotho.get(wrapper);

clotho.get('UCB');

//Query for objects matching args, then grab the first one
var listy = clotho.query({"city" : "Townsville"});
var existing = listy[0];
clotho.say("Clotho found " + existing.name);

//Create a partial replacement object
var args = {};
args.id = existing.id;
args.city = "Paris";

//Call set on the Object and print out its returned and persisted value
var result = clotho.set(args);
clotho.say('result is now: ' + JSON.stringify(result));

//Modify a complete existing object and set that as the new value
existing.name = "Party School"; 
existing.city = "New paris"; 
existing.state = "CA"; 
result = clotho.set(existing);
clotho.say('result is now: ' + JSON.stringify(result));

//Remove the ID field and create a new object
existing.id = undefined;
existing.name = "Creation College";
existing.city = "Clothoton";
existing.state = "Maine";
result = clotho.create(existing);
clotho.say('result is now: ' + JSON.stringify(result));


Test of say operations
clotho.say('hi there');

var obj = {};
obj.value = 'value is an acceptable token';
clotho.say(obj);

var obj = {};
obj.msg = 'msg, value, message are all acceptable tokens';
clotho.say(obj);

var obj = {};
obj.message = 'msg, value, message, args are all acceptable tokens';
clotho.say(obj);

obj.message = 'can stuff severity in the object';
obj.severity = 'text-error';
clotho.say(obj);

var cats = {};
cats.name = 'Fluffy';
clotho.say(cats);

clotho.say('this one doesnt pick up the second argument right now', 'text-error');

//Testing out clear, last say should fail
var cats = 'cats';
clotho.say(cats);
clotho.clear();
clotho.say(cats);


//Create a bunch of institutions with state Iowa then delete them
var i=0;
while(i<5) {
   clotho.create( {"lastModified":{"$date":"2013-06-13T02:47:51.288Z"},"name":"Temporary College","state":"Iowa","className":"org.clothocad.model.Institution","isDeleted":false,"country":"United States of America","city":"Townsville"} );
   i++;
}

var results = clotho.query({"state":"Iowa"});

clotho.destroy(results);

clotho.query({"state":"Iowa"});

*/
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
