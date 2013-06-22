println('About to run MindTest');

//Quick housekeeping before tests
var totoss = clotho.query({"city" : "Baltizam"});
if(totoss.length>0) {
    clotho.destroy(totoss);
}
totoss = clotho.query({"city" : "Whaletown"});
if(totoss.length>0) {
    clotho.destroy(totoss);
}

println('mind test runninig');

var newobj = clotho.create( {"name":"UCM","state":"MA","className":"org.clothocad.model.Institution","isDeleted":false,"country":"United States of America","city":"Baltizam"} );

println('mind test create works');

//Get the newobj by a get call (with different sloppiness)
var result = clotho.get(newobj.id);
if(result.name != "UCM") {
    throw new Exception();
}

println('mind test basic get works');

result = clotho.get(newobj);
if(result.name != "UCM") {
    throw new Exception();
}

var wrapper = {};
wrapper.data = newobj;
result = clotho.get(wrapper);
if(result.name != "UCM") {
    throw new Exception();
}


wrapper = {};
wrapper.data = newobj.id;
result = clotho.get(wrapper);
if(result.name != "UCM") {
    throw new Exception();
}

wrapper = [];
wrapper[0] = newobj;
result = clotho.get(wrapper);
if(result.name != "UCM") {
    throw new Exception();
}

result = clotho.get('UCM');
if(result.name != "UCM") {
    throw new Exception();
}

println('mind test sloppy get works');

//Query for objects matching args, then grab the first one
var listy = clotho.query({"city" : "Baltizam"});
var existing = listy[0];
if(!existing) {
    throw new Exception();
}

println('mind test basic query works');

//Create a partial replacement object, then call set
var args = {};
args.id = existing.id;
args.city = "Paris";
var result = clotho.set(args);
if(result.city != "Paris") {
    throw new Exception();
}

println('mind test basic set works');

//Modify a complete existing object and set that as the new value
existing.name = "Shamoo University"; 
existing.city = "Whaletown"; 
existing.state = "NR"; 
result = clotho.set(existing);
if(result.city != "Whaletown") {
    throw new Exception();
}
if(result.state != "NR") {
    throw new Exception();
}

println('mind test basic set 2 works');

//Test of say operations
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

println('mind test say works');

//clotho.say('this one doesnt pick up the second argument right now', 'text-error');

//Test query and destroy (again) to cleanup
var finalSet = clotho.query({"city" : "Whaletown"});
println("finalSEt" + finalSet);
if(finalSet.length!=1) {
    throw new Exception();
}

clotho.destroy(finalSet);
finalSet = clotho.query({"city" : "Whaletown"});
if(finalSet.length!=0) {
    throw new Exception();
}