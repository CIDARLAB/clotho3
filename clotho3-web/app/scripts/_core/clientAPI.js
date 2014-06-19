/**
 * @name ClientAPI
 *
 * @description
 * This is the client Clotho API - commands issued BY the server to be run on the client
 */
angular.module('clotho.core').service('ClientAPI',
	function (PubSub, Collector, Debug, $q, $injector, $window, $templateCache, $http, $rootScope, $location, $document) {

		var Debugger = new Debug('ClientAPI', '#dd99dd');

		//check for presence so don't rely on module
		if ( $injector.has('$clothoModal') ) {
			var $clothoModal = $injector.get('$clothoModal');
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
		 * @param {object} uuid UUID of view to show
		 * @param targetSelector CSS Selector of target to append to, otherwise widget area at bottom

		 note CAVEATS:
		 - currently, controllers etc. must be tied to Application.Extensions.___
		 */
		var display = function clientAPIDisplay(uuid, targetSelector) {

			var showDiv = angular.element('<div clotho-show="' + uuid + '"></div>');
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

			 //NOTE - jquery reliance (not in use)
			 target = data.target && $($clotho.appRoot).has(data.target) ? data.target : $($clotho.appRoot).find('[ng-view]');

			 $rootScope.$safeApply($http.get(template, {cache: $templateCache})
			 .success(function(precompiled) {

			 $clotho.extensions.mixin([dependencies, controller], angular.element(precompiled).appendTo(target), args)
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
			var el = angular.element(document.querySelector('[clotho-show="' + uuid + '"]'));

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
			Debugger.info(msg);
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

			//parse message text (before extend i.e. null -> undefined)
			if (angular.isString(data.text)) {
				//ok
			} else if (angular.isNumber(data.text)) {
				data.text = parseInt(data.text);
			} else if (data.text === null) {
				data.text = 'null';
			}	else if (angular.isUndefined(data.text)) {
				data.text = 'undefined';
			} else if (data.text === true) {
				data.text = 'true';
			} else if (data.text === false) {
				data.text = 'false';
			}

			var defaults = {
				'class': 'muted',
				'from': 'server',
				'timestamp': Date.now().valueOf()
			};

			data = angular.extend({}, defaults, data);

			//parse class
			var classMap = {
				success: 'success',
				warning: 'warning',
				error : 'danger',
				failure: 'danger',
				normal: 'success',
				muted: 'muted',
				info: 'info'
			};
			data.class = classMap[angular.lowercase(data.class)];

			PubSub.trigger('activityLog', data);
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


			angular.isDefined('$clothoModal') && $clothoModal.create({
				title : "Clotho Alert",
				content : "msg"
			})
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

			//todo - create clotho modal and show it

		};

		/**
		 * @name clientAPI.revisions
		 *
		 * @param {string} uuid
		 * @param {object} data
		 *
		 * @description
		 * Publish list of versions for a given resource on 'revisions:<uuid>'
		 */
		var revisions = function clientAPIRevisions(uuid, data) {
			PubSub.trigger('revisions:' + uuid, data);
		};

		/**
		 * @name clientAPI.changeUrl
		 *
		 * @param newUrl URL to change to, angular version
		 */
		var changeUrl = function (newUrl) {
			$location.path(newUrl);
		};

		/**
		 * @name clientAPI.edit
		 * @param uuid UUID of sharable to edit, opens in modal
		 */
		var edit = function (uuid) {
			if ($injector.has('$route')) {
				$location.path("/editor?id=" + uuid);
			} else {
				//need to go to full instance, this just has API or Command build
				say({
					text : 'Editor not available',
					from: 'client',
					class : 'error'
				});
			}
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
			$location.path('/trails/' + uuid);
		};

		return {
			collect: collect,
			broadcast: broadcast,
			log: log,
			say: say,
			alert: alert,
			display: display,
			hide: hide,
			help: help,
			revisions: revisions,

			//using routes
			changeUrl: changeUrl,
			edit: edit,
			startTrail: startTrail
		}
	});