/* Widget API */

var collector = new Object();

collector.add = function(data) {
    collector[data.id] = data;
}

collector.get = function(uuid) {
    return collector[uuid];
}