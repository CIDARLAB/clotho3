package org.clothocad.core.layers.communication;

public enum Channel {

    autocomplete, //Return potential commands that start with this substring
    autocompleteDetail, //Return the metadata for a Sharable
    submit, //Run this sloppy or concrete command
    login, //Log me into Clotho on this client with this login/password
    logout, //Log me out of Clotho
    changePassword, //Change my password to this new value

    log, //write a message to the server log for the Doo managing the current process
    say, //Post a message to all search bars for this Person
    
    get, //Get the Sharable with this uuid
    set, //Change the state of the Sharable with this uuid with this new JSON data
    query, //Find all Sharables that satisfy this formally-expressed constraint
    create, //Create a new sharable with this Schema

    run, //run this Function
    learn, //associate this String with this execution statement
    show, //Instantiate the View with this uuid
    startTrail, //Show this trail
    destroy, //Delete a sharable with this uuid
    edit, //Open a sharable in an editor
    listen, //Listen to a pubsub channel for events, and in response do this callback
    unlisten, //Remove a listener for an event
    alert, //Post an alert to this Person, wherever they say they want to be contacted in preferences.  This isn't a javascript alert().  That is just one implementation.


}
