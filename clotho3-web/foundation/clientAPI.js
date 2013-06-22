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
 */
Application.Foundation.service('ClientAPI', ['PubSub', 'Collector', '$q', '$templateCache', '$http', '$rootScope', '$location', '$compile', '$dialog', function(PubSub, Collector, $q, $templateCache, $http, $rootScope, $location, $compile, $dialog) {

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
     * @name clientAPI.collect
     *
     * @param {object} data
     * JSON as outlined below:
     *
     *  {
     *      "uuid" : {string} UUID
     *      "type" : {string} json | html | css | javascript
     *      "model" : {JSON | string} for json, object to store. for url, url as string
     *      "isURL" : if URL, true if URL, false if template
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

    var edit = function(uuid) {
        var dialog_opts = {
            backdrop: true,
            keyboard: true,
            backdropClick: true,
            template:  '<form sharable-editor name="sharableEditor" uuid="'+uuid+'" class="span6 form-horizontal well" novalidate></form>'
        };
        var d = $dialog.dialog(dialog_opts);
        d.open();
    };

    var changeUrl = function(newUrl) {
        $location.path(newUrl);
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
    var display_old = function clientAPIDisplay(uuid, args) {

        //as a simple demo, let's assume this is all that is sent
        var model = args.model;
        var template = args.template;

    };

    /**
     * @name clientAPI.display
     *
     * @param {object} data
     * format:
        {
            "template" : <url>,         // required
            "target" : <DOM ELEMENT>    // suggested, or absolute positioning in CSS
            "args" : {<object>}         // data to copy onto $scope
            "controller" : <url>,       // optional
            "dependencies" : [
                <urls>                  // required if in controller
            ],
            styles : {
                <styles>
                //e.g.
                "background-color" : "#FF0000"
            }
        }

     note CAVEATS:
     - currently, controllers etc. must be tied to Application.Extensions.___
     */
    var display = function clientAPIDisplaySimple(data) {

        //console.log(data);

        var template = data.template,
            controller = data.controller || "",
            args = data.args || {},
            dependencies = data.dependencies || [],
            styles = data.styles || {},
            target = data.target && $($clotho.appRoot).has(data.target) ? data.target : $clotho.appRoot;

        $rootScope.$safeApply($http.get(template, {cache: $templateCache})
            .success(function(precompiled) {

                Application.mixin([dependencies, controller], $(precompiled).appendTo(target), args)
                    .then(function(div) {
                        //testing
                        //console.log($position.position(div));
                        //console.log($position.position(target));
                        div.css(styles);
                    });
            })
            .error(function (result) {
                console.log("error getting template");
            })
        );
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
        console.log('Hit say');
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

        PubSub.trigger('serverAlert', {});
        $rootScope.$safeApply($dialog.serverAlert(msg)
            .open()
            .then(function(result){
                console.log('dialog closed with result: ' + result);
            })
        );
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
     * @param {object} data
     *
     * @description
     * Publish list of versions for a given resource on "revisions:<uuid>"
     */
    var revisions = function clientAPIRevisions(uuid, data) {
        PubSub.trigger('revisions:'+uuid, data);
    };

    /**
     * @name Clotho.startTrail
     *
     * @param {string} uuid
     *
     * @description
     * start a trail with a given uuid
     */
    var startTrail = function clothoAPI_startTrail(uuid) {
        $location.path("/trails/" + uuid);
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
        console.log('Hit autocomplete');
        console.log(JSON.stringify(list));

        PubSub.trigger('autocomplete', list);
    };

    var autocompleteDetail = function clientAPIAutocompleteDetail(obj) {

        console.log('Hit autocompletedetail');
        console.log(obj);
        
        var id = obj.command_object.function_id;

        Collector.storeModel("detail_" + id, obj);
        PubSub.trigger('autocompleteDetail_'+id, obj);
    };


    return {
        mapCommand : mapCommand,
        collect : collect,
        edit : edit,
        changeUrl : changeUrl,
        notify : notify,
        broadcast : broadcast,
        log : log,
        say : say,
        alert : alert,
        display : display,
        display_simple : display,
        display_old : display_old,
        hide : hide,
        help : help,
        revisions : revisions,
        startTrail : startTrail,
        autocomplete : autocomplete,
        autocompleteDetail: autocompleteDetail
    }
}]);