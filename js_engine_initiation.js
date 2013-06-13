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

*/
var clotho = {};

//var complete = clotho.autocomplete("walk the");
//clotho.say(JSON.stringify(complete));
clotho.autocomplete = function(args, callback) {
    var sresult = clothoJava.autocomplete(args);
    var result = JSON.parse(sresult);
    if(callback) callback(result );
    else return result;
};

clotho.autocompleteDetail = function(args, callback) {
    var sresult = clothoJava.autocompleteDetail(args);
    var result = JSON.parse(sresult);
    if(callback) callback(result );
    else return result;
};

clotho.submit = function(args, callback) {
    clothoJava.submit(args);
    if(callback) callback( );
};

clotho.learn = function(nativeCmd, jsCmd, callback) {
    clothoJava.learn(nativeCmd, jsCmd);
    if(callback) callback( );
};

//clotho.say("Says come up in the activity log");
clotho.say = function(args, callback) {
    clothoJava.say(args);
    if(callback) callback();
};

clotho.alert = function(args, callback) {
    clothoJava.alert(args);
    if(callback) callback();
};

clotho.log = function(args, callback) {
    clothoJava.log(args);
    if(callback) callback();
};

clotho.note = function(args, callback) {
    clothoJava.note(args);
    if(callback) callback();
};

//var paris = clotho.get('51b8b90050765ba45b5a8217');
//clotho.say(paris.state);
clotho.get = function(args, callback) {
    var sresult = clothoJava.get(args);

    var result = JSON.parse(sresult);

    if(callback) callback(result);
    else return result;
};




clotho.set = function(value, callback) {
    var jsonstr = JSON.stringify(value);
    var sresult = clothoJava.set(jsonstr);

    var result = JSON.parse(sresult);

    if(callback) callback(result);
    else return result;
};

//var inst = {"city":"Paris","className":"org.clothocad.model.Institution","country":"United States of America","isDeleted":false,"lastModified":{"$date":"2013-06-08T02:48:10.254Z"},"name":"Test institution","state":"Massachusetts"};clotho.say(inst.city);   var result = clotho.create(inst);
clotho.create = function(args, callback) {
    var jsonstr = JSON.stringify(args);
    var sresult = clothoJava.create(jsonstr);
    var result = JSON.parse(sresult);
    if(callback) callback(result );
    else return result;
};

clotho.destroy = function(sharableId) {
    clothoJava.destroy(sharableId);
    if(callback) callback( );
};

//clotho.query('city', 'Townsville');
clotho.query = function(query, callback) {
    var jsonstr = JSON.stringify(query);
    var slisty = clothoJava.query(jsonstr);
    //println('query has slisty as ' + slisty);
    var listy = JSON.parse(slisty);
    if(callback) callback(listy);
    else return listy;
};

