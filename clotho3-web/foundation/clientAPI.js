'use strict';

/**
 * @name clientAPI
 *
 * @description
 * This is the client Clotho API - commands issued BY the server to be run on the client
 *
 * Notes:
 *  - not really sure at this point which functions should `return` and which should have a `callback` ... so some of the functions have a parameter callback, particularly those which are asynchronous, in the event we want to pass data back, or allow for more custom operations to be made using the data passed to / from the command.
 *
 *  FUTURE - once the Command Bar has been integrated with Angular, API commands to interact with it
 *
 */
Application.Foundation.service('ClientAPI', ['PubSub', 'Collector', '$templateCache', '$http', '$rootScope', '$location', function(PubSub, Collector, $templateCache, $http, $rootScope, $location) {

    /**
     * @name clientAPI.mapCommand
     *
     * @param command {string} the command to run
     * @param data {object} JSON to pass along
     *
     * @description
     * Given a command
     *
     */
    var mapCommand = function clientAPIMapCommand(command, data) {
        //todo - should have more logic than simply calling ClientAPI[command](data)
        console.log("CLIENTAPI\tmapCommand called on channel " + command);

    };


    /**
     * todo : think this belongs in ServerAPI
     * @name clientAPI.clone
     *
     * @param {string} uuid
     *
     * @description
     * copy a sharable, which won't be updated later
     *
     * @return
     * returns a full JSON representation of a sharable
     */

    var clone = function clientAPIClone(uuid) {

    };

    /**
     * @name clientAPI.collect
     *
     * @param {object} data
     * JSON as outlined below:
     *
     *  {
     *      "uuid" : {string} UUID
     *      "type" : {string} json | html | css | javascript
     *      "data" : {JSON} for json, object to store
     *      "url" : {string} for non-json, resource to collect
     *  }
     *
     * @description
     * Store a resource
     * - Models (JSON) will go into the collector
     * - templates (HTML) will be cached in the Angular template cache
     * - CSS, JS, will be downloaded and cached automatically with the right headers
     *      - they should be versioned e.g. clotho.com/css/some_css.css?version=1.1
     *
     * Because models will be updated automatically, collecting a model will update a view
     *
     */
    //future - share templatecache across apps
    //todo - better handling of url vs. model for template -- use regex and get rid of isURL

    var collect = function clientAPICollect(data) {

        var type = data.type;
        var uuid = data.uuid;
        var model_or_url = data.model;
        var isURL = data.isURL;

        switch (type) {
            case "json" :
                Collector.storeModel(uuid, model_or_url);
                break;
            case "html" :
                if (isURL) {

                    //todo - remove and reset if exists
                    // http://deansofer.com/posts/view/14/AngularJs-Tips-and-Tricks-UPDATED#refresh-partials
                    /*
                     //check cache
                     var template = $templateCache.get(model_or_url);
                     console.log(template);
                     if (template) {
                     //don't do anything
                     console.log("already exists");
                     } else {
                     $http.get(model_or_url, {cache:$templateCache})
                     .then(function(result) {
                     console.log("got template: " + uuid);
                     console.log(result);
                     })
                     }*/


                    $rootScope.$safeApply($http.get(model_or_url, {cache:$templateCache}));

                } else {
                    //sending the actual template (isURL = false)
                    $templateCache.put(uuid, model_or_url);
                }
                break;
            case "js-template" :
                //todo - handle passing template as script
                break;
            case "js" :
                $script(model_or_url, uuid);
                break;
            default :
                console.log("CLIENTAPI\tcollect\ttype not sent just download: " + uuid);
                $rootScope.$safeApply($http.get(model_or_url));
                break;
        }
    };

    //note - temporary - need to integrate modal
    var edit = function(uuid) {
        $location.path('/editor/' + uuid)
    };

    /**
     * @name clientAPI.notify
     *
     * @param {string} uuid UUID of model to update
     * @param {object} obj  Model
     * @param {function} callback
     *
     * @description
     * ??
     * concerns -- I think it's best to pass in a callback in serverAPI.run() (and every other method) so you can pass in $scope etc. without having to dig through the angular application looking for it, or angular applications on a given page. Communication between apps should happen on a separate channel.
     *
     */
    var notify = function clientAPINotify(uuid, obj, callback) {

    };

    /**
     * @name clientAPI.revert
     *
     * @param {string} uuid UUID of resource with unwanted changes
     * @param {string} timestamp Date of desired version
     * @param {function} callback
     *
     * @description
     * Revert a model to an earlier version
     *
     * @returns {object} version at timestamp of resource with passed UUID
     *
     */
    var revert = function clientAPIRevert(uuid, timestamp, callback) {

    };


    /**
     * @name clientAPI.broadcast
     *
     * @param {object} obj  Object to pass to PubSub
     *  Of the form:
     *  {
     *      channel: <channel>,
     *      data: <data>
     *  }
     *
     * @description
     * Publish an event to PubSub
     *
     */
    var broadcast = function clientAPIBroadcast(obj) {
       PubSub.trigger(obj.channel, obj.data);
    };

    /**
     * @name clientAPI.display
     *
     * @param {string} uuid
     * @param {object} args
     *  - ? [custom]
     *  - position
     *      - parent : { uuid : [uuid of div to insert into] }
     *      - absposition : { {int} x, {int} y }
     *      - mode : { mode }

     *
     * @description
     * Show a view on the client
     *
     */
    var display = function clientAPIDisplay(uuid, args) {

        //as a simple demo, let's assume this is all that is sent
        var model = args.model;
        var template = args.template;

    };

    /**
     * @name clientAPI.hide
     *
     * @param {string} uuid
     * @param {function} callback
     *
     * @description
     * Hide a view on the client
     *
     */
    var hide = function clientAPIHide(uuid, callback) {

    };

    /**
     * @name clientAPI.log
     *
     * @param {string} msg
     *
     * @description
     * Writes a message to the console
     *
     */
    var log = function clientAPILog(msg) {
        console.log(msg);
    };

    /**
     * @name clientAPI.say
     *
     * @param {object} data
     * {"msg" : <msg>}
     *
     * @description
     * Adds a message to the Command Bar's Activity Log
     *
     */
    var say = function clientAPISay(data) {
        PubSub.trigger("activityLog", data);
    };

    /**
     * @name clientAPI.alert
     *
     * @param {string} msg
     *
     * @description
     * Alerts a message
     *
     */
    var alert = function clientAPIAlert(msg) {
        console.log(msg);
        window.alert(msg);
    };

    /**
     * @name clientAPI.help
     *
     * @param {string} uuid .... what parameter should be sent?
     *
     * @description
     * Get help for a mode or instance or something
     */
    var help = function clientAPIHelp(uuid) {

    };

    /**
     * @name clientAPI.revisions
     *
     * @param {string} uuid
     *
     * @description
     * Get a list of versions for a given resource
     *
     * @returns {array} Array containing all (?) timestamps as a list of strings
     */
    var revisions = function clientAPIRevisions(uuid) {

    };

    // ---- COMMAND BAR ----

    /**
     * @name clientAPI.autocomplete
     *
     * @param {array} list Array of autocompletions
     *
     * @description
     * Publishes autocompletions to PubSub for listeners to pick up
     */
    var autocomplete = function clientAPIAutocomplete(list) {
        PubSub.trigger('autocomplete', list);
    };

    /*var autocompleteDetail = function clientAPIAutocompleteDetail(obj) {
        Collector.storeModel("detail_" + obj.uuid, obj);
    };*/


    return {
        mapCommand : mapCommand,
        clone : clone,
        collect : collect,
        edit : edit,
        notify : notify,
        revert : revert,
        broadcast : broadcast,
        log : log,
        say : say,
        alert : alert,
        display : display,
        hide : hide,
        help : help,
        revisions : revisions,
        autocomplete : autocomplete,
        //autocompleteDetail: autocompleteDetail
    }
}]);