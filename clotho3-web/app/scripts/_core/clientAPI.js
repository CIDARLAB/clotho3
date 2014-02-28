/**
 * @name ClientAPI
 *
 * @description
 * This is the client Clotho API - commands issued BY the server to be run on the client
 */
angular.module('clotho.core').service('ClientAPI',
	function(PubSub, Collector, $q, $templateCache, $http, $rootScope, $location, $compile) {

		//todo - verify this works and retrieves the service properly
		var interfaceModulePresent = false,
			$modal;
			try {
				angular.module('clotho.interface');
				interfaceModulePresent = true;
				$modal = angular.injector('clotho.interface').get('$modal');
			} catch (err) {
				// not present
			}

    /**
     * @name clientAPI.collect
     *
     * @param {object} data with field id
     *
     * @description
     * Store an object. Generally, objects are requested via get(). However, the server can add objects to the client collector via collect().
     */

    var collect = function clientAPICollect(data) {
        var model = data,
            id = model.id || model.uuid || false;

        Collector.storeModel(id, model);
    };

    /**
     * @name clientAPI.edit
     * @param uuid UUID of sharable to edit, opens in modal
     */
    var edit = function(uuid) {
	    if (interfaceModulePresent) {
		    var dialog_opts = {
			    backdrop: true,
			    keyboard: true,
			    backdropClick: true,
			    templateUrl:  '<form sharable-editor ng-model="'+uuid+'" class="col-sm-6 form-horizontal well"></form>'
		    };
		    $modal.open(dialog_opts);
	    } else {
		    $location.path('/editor/' + uuid)
	    }
    };

    /**
     * @name clientAPI.changeUrl
     *
     * @param newUrl URL to change to, angular version
     */
    var changeUrl = function(newUrl) {
        $location.path(newUrl);
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
       PubSub.trigger(obj.channel, [obj.data]);
    };

    /**
     * @name clientAPI.display_old
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
     * @param {object} uuid UUID of view to show
     * @param targetSelector CSS Selector of target to append to, otherwise widget area at bottom

     note CAVEATS:
     - currently, controllers etc. must be tied to Application.Extensions.___
     */
    var display = function clientAPIDisplay(uuid, targetSelector) {

	    var showDiv = angular.element('<div clotho-show="'+uuid+'"></div>');
	    var targetEl = document.querySelector(targetSelector);

	    if (!targetEl) {
		    targetEl = document.getElementById('clothoAppWidgets');
	    }

	    targetEl = angular.element(targetEl);

	    targetEl.append($compile(showDiv)($rootScope));

	    /*
	    console.log(data);

        var template = data.template,
            controller = data.controller || "",
            args = data.args || {},
            dependencies = data.dependencies || [],
            styles = data.styles || {},
            target = data.target && $($clotho.appRoot).has(data.target) ? data.target : $($clotho.appRoot).find('[ng-view]');

        $rootScope.$safeApply($http.get(template, {cache: $templateCache})
            .success(function(precompiled) {

                $clotho.extensions.mixin([dependencies, controller], $(precompiled).appendTo(target), args)
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
        */
    };

    /**
     * @name clientAPI.hide
     *
     * @param {string} uuid of view to remove
     * @param {function} callback which passed the removed element
     *
     * @description
     * Hide a view on the client
     *
     */
    var hide = function clientAPIHide(uuid, callback) {
	    var el = angular.element(document.querySelector('[clotho-show="'+uuid+'"]'));

	    callback.apply(null, el.remove());
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
     * {
            "text" : msg,
            "from" : sender,
            "class" : css,
            "timestamp" : timestamp
       }
     *
     * @description
     * Adds a message to the Command Bar's Activity Log
     *
     */
    var say = function clientAPISay(data) {
	      var defaults = {
		      'class' : 'success',
		      'from' : 'server',
		      'timestamp' : Date.now().valueOf()
	      };

	      angular.extend(defaults, data);

        PubSub.trigger("activityLog", [defaults]);
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

        PubSub.trigger('serverAlert');

		    if (interfaceModulePresent) {
			    $rootScope.$safeApply($modal.serverAlert(msg)
				    .result
				    .then(function(result){
					    console.log('dialog closed with result: ' + result);
				    })
			    );
		    } else {
			    window.alert(msg);
		    }
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
        PubSub.trigger('revisions:'+uuid, [data]);
    };

    /**
     * @name clientAPI.startTrail
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

        PubSub.trigger('autocomplete', [list]);
    };

    var autocompleteDetail = function clientAPIAutocompleteDetail(obj) {

        console.log('Hit autocompletedetail');
        console.log(obj);
        
        var id = obj.command_object.function_id;

        Collector.storeModel("detail_" + id, obj);
        PubSub.trigger('autocompleteDetail_'+id, [obj]);
    };


    return {
        collect : collect,
        edit : edit,
        changeUrl : changeUrl,
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
        //autocomplete : autocomplete,
        autocompleteDetail: autocompleteDetail
    }
});