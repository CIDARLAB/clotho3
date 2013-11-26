'use strict';

//future - keep list of components already mixed in

window.$clotho.extensions = $clotho.extensions = {};

angular.module('clotho.extensions', [])
.config(function($routeProvider, $controllerProvider, $compileProvider, $filterProvider, $provide) {

	$clotho.extensions.providers = {
		$routeProvider: $routeProvider,
		$controllerProvider: $controllerProvider,
		$compileProvider: $compileProvider,
		$filterProvider: $filterProvider,
		$provide: $provide
	};

	$clotho.extensions.getQueue = function() {
		return angular.module('clotho.extensions')._invokeQueue;
	};

	$clotho.extensions.registeredQueue = $clotho.extensions.getQueue().length;

})
.run(function($rootScope, $q, $timeout, $templateCache, $http, $rootElement) {

	//need to call this before compiling new element
	$clotho.extensions.processQueue = function() {
		var queue = $clotho.extensions.getQueue();

		for(var i=$clotho.extensions.registeredQueue;i<queue.length;i++) {
			var call = queue[i];
			// call is in the form [providerName, providerFunc, providerArguments]
			var provider = $clotho.extensions.providers[call[0]];
			if(provider) {
				// e.g. $controllerProvider.register("Ctrl", function() { ... })
				//console.log("provider exists");
				provider[call[1]].apply(provider, call[2]);
			}
		}

		$clotho.extensions.registeredQueue = $clotho.extensions.getQueue().length;
	};

	$clotho.extensions.recompile = function(element, args) {
		//can't compile already-compiled elements or cause problems
		if (typeof element == 'undefined') {return;}
		args = args || {};

		//todo - check for class ng-scope in what compile -- don't wanna recompile

		//todo - should inherit from closest relative scope - not rootscope

		$rootElement.injector().invoke(function($compile, $rootScope) {
			var scope = $rootScope.$new();
			angular.extend(scope, args);
			$compile($(element))(scope);
			$rootScope.$apply();
		});
	};

	/**
	 * @name Application.mixin
	 *
	 * @param urls URLs of dependencies. Only downloaded if hasn't been already
	 * @param element Element to be compiled. Necessary for compiling (to avoid recompiling)
	 * @param args Arguments to extend the scope
	 * @returns {promise} Element Promise to be fulfilled on successful addition, value is element
	 */
	$clotho.extensions.mixin = function(urls, element, args) {

		if (angular.isUndefined(urls) || urls == '') {
			//console.log('[Application.mixin] no url - return empty resolved promise');
			return $q.when('no mixin url');
		}

		var deferred = $q.defer();

		$script(urls, function() {
			$clotho.extensions.processQueue();
			$clotho.extensions.recompile(element, args);
			$rootScope.$safeApply(deferred.resolve(element));
		});

		return deferred.promise;
	};

	/**
	 * @name Application.script
	 *
	 * @description Downloads and executes a script, using cache-busting. Timestamp is appended to the script, to ensure it will run each time.
	 *
	 * @param {string|Array} urls URLs of scripts to be executed
	 * @returns {Promise} Promise which is resolved when all scripts have been executed
	 */
	$clotho.extensions.script = function(urls) {

		if (angular.isUndefined(urls) || urls == '') {
			return $q.when('no script url');
		}

		var deferred = $q.defer(),
			downloads; //don't want to overwrite source urls with timestamp

		if (angular.isString(urls))
			downloads = urls + '?_=' + Date.now();
		else {
			downloads = [];
			angular.forEach(urls, function(url) {
				downloads.push(url + '?_=' + Date.now());
			});
		}

		$script(urls, function() {
			$rootScope.$safeApply(deferred.resolve());
		});

		return deferred.promise;
	};

	/**
	 * @name Application.css
	 *
	 * @description Downloads and appends a CSS file to the page head, so that it will be applied properly
	 *
	 * @param {URL} url URL of CSS file
	 *
	 * @returns {Promise} Promise to be fulfilled once CSS files are downloaded and appended
	 */
	$clotho.extensions.css = function(url) {

		if (angular.isUndefined(url) || url == '') {
			return $q.when('no css url');
		}

		//todo - track so only added once
		//todo - handle array

		var deferred = $q.defer();

		if (document.createStyleSheet) { //IE
			document.createStyleSheet(url);
			$rootScope.$safeApply(deferred.resolve());
		} else {
			var link = document.createElement("link");
			link.type = "text/css";
			link.rel = "stylesheet";
			link.href = url;
			document.getElementsByTagName("head")[0].appendChild(link);
			$rootScope.$safeApply(deferred.resolve());
		}

		return deferred.promise;
	};

	/**
	 * @name Application.cache
	 *
	 * @description Downloads caches an angular template for later use
	 *
	 * @param {URL} url URL of angular template file
	 *
	 * @returns {Promise} Promise to be fulfilled once CSS files are downloaded and appended
	 */
	$clotho.extensions.cache = function(url) {

		if (angular.isUndefined(url) || url == '') {
			return $q.when();
		}

		var deferred = $q.defer();

		$http.get(url, {cache:$templateCache})
			.success(function() {deferred.resolve()})
			.error(function() {deferred.reject()});

		return deferred.promise;
	};

	/**
	 * @name Application.bootstrap
	 * @previous Clotho.bootstrap
	 *
	 * @param appInfo {object} Object with necessary information to bootstrap, minimally including:
	 * {
     *      "moduleName" : <Name of module as defined in Angular>
     *      "moduleUrl" : "<URL to module js file>",
     * }
	 *
	 * @returns {Promise} Array of selectors in form: [<appUUID>, <jQuery Selector>]
	 * @description
	 * Load a widget and bootstrap it. appInfo must contain a full module. for simply adding components to the stack, use mixin()
	 */
	var widgetID = 0;
	$clotho.extensions.bootstrap = function (appInfo) {
		widgetID++;

		var deferred = $q.defer();

		//angular version
		//note angular returns parent, not appended element
		//todo - if want this, select appropriate child element
		//var insertInto = angular.element(document).find("ng-app-clothoWidgets").append(angular.element('<div clotho-widget clotho-widget-uuid="'+appUUID+'" clotho-widget-name="'+appInfo.moduleName+'"></div>').append('<div ng-view></div>'));


		//todo - move away from ng-view and routing - just use a template


		//jQuery version
		var insertInto = $($('<div clotho-widget clotho-widget-uuid="'+widgetID+'" clotho-widget-name="'+appInfo.moduleName+'"></div>').append('<div ng-view></div>')).appendTo($clotho.appWidgets);

		$clotho.extensions.script(appInfo.moduleUrl).then(function() {
			angular.bootstrap(insertInto, [appInfo.moduleName]);
			deferred.resolve([widgetID, "[clotho-widget-uuid="+widgetID+"]"]);
		});

		return deferred.promise;
	}

});
