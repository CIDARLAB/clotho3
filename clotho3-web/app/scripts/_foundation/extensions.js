/*!
 Async script loader
 $script.js v1.3
 usage: https://github.com/ded/script.js
 */
!function(a,b,c){function t(a,c){var e=b.createElement("script"),f=j;e.onload=e.onerror=e[o]=function(){e[m]&&!/^c|loade/.test(e[m])||f||(e.onload=e[o]=null,f=1,c())},e.async=1,e.src=a,d.insertBefore(e,d.firstChild)}function q(a,b){p(a,function(a){return!b(a)})}var d=b.getElementsByTagName("head")[0],e={},f={},g={},h={},i="string",j=!1,k="push",l="DOMContentLoaded",m="readyState",n="addEventListener",o="onreadystatechange",p=function(a,b){for(var c=0,d=a.length;c<d;++c)if(!b(a[c]))return j;return 1};!b[m]&&b[n]&&(b[n](l,function r(){b.removeEventListener(l,r,j),b[m]="complete"},j),b[m]="loading");var s=function(a,b,d){function o(){if(!--m){e[l]=1,j&&j();for(var a in g)p(a.split("|"),n)&&!q(g[a],n)&&(g[a]=[])}}function n(a){return a.call?a():e[a]}a=a[k]?a:[a];var i=b&&b.call,j=i?b:d,l=i?a.join(""):b,m=a.length;c(function(){q(a,function(a){h[a]?(l&&(f[l]=1),o()):(h[a]=1,l&&(f[l]=1),t(s.path?s.path+a+".js":a,o))})},0);return s};s.get=t,s.ready=function(a,b,c){a=a[k]?a:[a];var d=[];!q(a,function(a){e[a]||d[k](a)})&&p(a,function(a){return e[a]})?b():!function(a){g[a]=g[a]||[],g[a][k](b),c&&c(d)}(a.join("|"));return s};var u=a.$script;s.noConflict=function(){a.$script=u;return this},typeof module!="undefined"&&module.exports?module.exports=s:a.$script=s}(this,document,setTimeout);

//todo - add routes dynamically (forces ngRoute dependency)

//note - this module should not contain any components (controller, directives, etc.) to be loaded at initial bootstrap. $inject calls are re-routed through provider to add them manually

//note - you can always add a new module, and bootstrap a new widget. That widget can inherit the functionality of the clotho webapp by including the appropriate package (at time of writing, clotho.fullPackage). Then you can have config and run blocks, ignore routing, run autonomously, configure providers, etc.

angular.module('clotho.extensions', [])
.config(function($controllerProvider, $compileProvider, $filterProvider, $provide) {

	window.$clotho.extensions = $clotho.extensions = {};

	$clotho.extensions.providers = {
		$controllerProvider: $controllerProvider,
		$compileProvider: $compileProvider,
		$filterProvider: $filterProvider,
		$provide: $provide
	};

	// adapted from Ifeanyi Isitor: http://ify.io/lazy-loading-in-angularjs/
	// preserve the original references. angular.module('clotho.extensions').controller() is the same as $clotho.extensions._controller(). You can use these references and process the queue manually. Recommended you use the provider-based interface below though, and download components which take the form $clotho.extensions.controller() [without the underscore] which will be added automatically
	$clotho.extensions._controller = $clotho.extensions.controller;
	$clotho.extensions._service = $clotho.extensions.service;
	$clotho.extensions._factory = $clotho.extensions.factory;
	$clotho.extensions._value = $clotho.extensions.value;
	$clotho.extensions._directive = $clotho.extensions.directive;


	// Convert controllers, services, factories, directives, filters to provider-based alternatives, so that when they are loaded later, they are automatically registered
	$clotho.extensions.controller = function( name, constructor ) {
		$controllerProvider.register( name, constructor );
		return( this );
	};
	$clotho.extensions.service = function( name, constructor ) {
		$provide.service( name, constructor );
		return( this );
	};
	$clotho.extensions.factory = function( name, factory ) {
		$provide.factory( name, factory );
		return( this );
	};
	$clotho.extensions.value = function( name, value ) {
		$provide.value( name, value );
		return( this );
	};
	$clotho.extensions.directive = function( name, factory ) {
		$compileProvider.directive( name, factory );
		return( this );
	};
	$clotho.extensions.filter = function( name, factory ) {
		$filterProvider.filter( name, factory );
		return( this );
	};

})
.run(function($rootScope, $q, $timeout, $templateCache, $http, $rootElement, $compile) {

	//these functions are for mixing in components that are not added via the provider interface above
	var getQueue = function() {
		return angular.module('clotho.extensions')._invokeQueue;
	};
	var registeredQueue = getQueue().length;
	//Only needs to be called when not using the provider interface above, i.e. components that have been prefaced, e.g. angular.module('clotho.extensions').controller() instead of $clotho.extensions.controller()
		//todo - verify that processing queue after adding components using provider-based mechanism doesn't cause problems with queue length mismatch
	var processQueue = function() {
		var queue = getQueue();

		for(var i=registeredQueue;i<queue.length;i++) {
			var call = queue[i];
			// call is in the form [providerName, providerFunc, providerArguments]
			var provider = $clotho.extensions.providers[call[0]];
			if(provider) {
				// e.g. $controllerProvider.register("Ctrl", function() { ... })
				//console.log("provider exists");
				provider[call[1]].apply(provider, call[2]);
			}
		}

		registeredQueue = getQueue().length;
	};





	//note - create a new isolate scope, deleting scope and children if present
	//todo - DEPRECATE. can't really scrub DOM of compiled code, so should use bootstrap for something that needs to be compiled (i.e. cache template and bootstrap that div, as easy to recompile a template in $templateCache)
	$clotho.extensions.recompile = function(element, args) {
		if (typeof element == 'undefined') {return;}
		args = args || {};

		if (element.hasClass('ng-scope')) element.scope().$destroy();

		$rootElement.injector().invoke(function($compile, $rootScope) {
			var scope = $rootScope.$new(true);
			angular.extend(scope, args);
			$compile($(element))(scope);
			$rootScope.$apply();
		});
	};




	$clotho.extensions.extend = angular.extend;

	//not really recommended for use...
	$clotho.extensions.extendPrimaryRootscope = function(args) {
		$clotho.extensions.extend($rootScope, args);
	};





	/**
	 * @name Application.mixin
	 *
	 * @description Will download URLs only once
	 *
	 * @param urls URLs of dependencies. Only downloaded if hasn't been already
	 * @returns {Promise} Promise to be fulfilled on successful addition, value is urls passed
	 */
	$clotho.extensions.mixin = function(urls) {

		if (angular.isUndefined(urls) || urls == '') {
			return $q.when('no mixin url');
		}

		var deferred = $q.defer();

		//timeout our requests at 5 seconds
		var timeoutPromise = $timeout(function() { deferred.reject(null) }, 5000);

		$script(urls, function() {
			$timeout.cancel(timeoutPromise);
			$rootScope.$safeApply(deferred.resolve(urls));
		});

		return deferred.promise;
	};

	/**
	 * @name Application.script
	 *
	 * @description Downloads and executes a script or scripts, using cache-busting. Timestamp is appended to the script, to ensure it will run each time.
	 *
	 * @param {string|Array} urls URLs of scripts to be executed
	 * @returns {Promise} Promise which is resolved when all scripts have been executed
	 */
	$clotho.extensions.script = function(urls) {

		if (angular.isUndefined(urls) || urls.length == 0) {
			return $q.when('no script url');
		}

		var downloads = []; //don't want to overwrite source urls with timestamp

		if (angular.isString(urls)) {
			urls = [urls];
		}

		angular.forEach(urls, function(url) {
			downloads.push(url + '?_=' + Date.now());
		});

		return $clotho.extensions.mixin(downloads);
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
	var css_urls_downloaded = [];
	$clotho.extensions.css = function(url) {

		if (angular.isUndefined(url) || url == '') {
			return $q.when('no css url');
		}

		if (_.indexOf(css_urls_downloaded, url) > -1)
			return $q.when('CSS url already added');

		//todo - handle array

		var deferred = $q.defer();

		//timeout our requests at 5 seconds
		var timeoutPromise = $timeout(function() { deferred.reject(null) }, 5000);

		if (document.createStyleSheet) { //IE
			document.createStyleSheet(url);
			$rootScope.$safeApply(deferred.resolve());
		} else {
			var link = document.createElement("link");
			link.type = "text/css";
			link.rel = "stylesheet";
			link.href = url;
			document.getElementsByTagName("head")[0].appendChild(link);
			$timeout.cancel(timeoutPromise);
			$rootScope.$safeApply(deferred.resolve());
		}
		css_urls_downloaded.push(url);

		return deferred.promise;
	};

	/**
	 * @name Application.cache
	 *
	 * @description Downloads caches an angular template for later use. Forces addition to the cache, under ID of the passed URL. So, to use later, e.g. use <div ng-include="url_you_passed"></div>. Note that the template is cached in the primary Clotho App, so to access it in a separately bootstrapped app, you'll need to list the appropriate angular module as a dependency
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

		//timeout our requests at 5 seconds
		var timeoutPromise = $timeout(function() { deferred.reject(null) }, 5000);

		$http.get(url)
			.success(function(data) {
				$timeout.cancel(timeoutPromise);
				$templateCache.put(url, data);
				deferred.resolve(data)
			})
			.error(function(data) {deferred.reject(data)});

		return deferred.promise;
	};


	/**
	 * @name $clotho.extensions.bootstrap
	 *
	 * @description This is just a reference to angular.bootstrap. Bootstraps a new app manually, creating a new (isolate) $rootScope outside the flow of the parent app. Need to define clotho angular modules for them to be present in this app if you want their functionality. Note that this method expects all dependencies to have already been downloaded
	 */
	$clotho.extensions.bootstrap = angular.bootstrap;








	$clotho.extensions.determineUrlExtension = function ( url ) {
		//The extension is always the last characters before the ? and after a period.
		//accounting for the possibility of a period in the query string
		var b = url.split('?')[0];
		return b.substr(b.lastIndexOf('.')+1);
	};



	//todo - incorporate separately? don't need to know type this way
	var headEl = document.getElementsByTagName('head')[0];

		/*
		EXAMPLE USAGE

		// note on angular site (docs app)
		 // dynamically add base tag as well as css and javascript files.
		 // we can't add css/js the usual way, because some browsers (FF) eagerly prefetch resources
		 // before the base attribute is added, causing 404 and terribly slow loading of the docs app.

		 var jquery = true;
		 addTag('base', {href: baseUrl});
		 addTag('link', {rel: 'stylesheet', href: 'components/bootstrap/css/' + (debug ? 'bootstrap.css' : 'bootstrap.min.css'), type: 'text/css'});
		 addTag('link', {rel: 'stylesheet', href: 'components/font-awesome/css/' + (debug ? 'font-awesome.css' : 'font-awesome.min.css'), type: 'text/css'});
		 addTag('link', {rel: 'stylesheet', href: 'css/prettify.css', type: 'text/css'});
		 addTag('link', {rel: 'stylesheet', href: 'css/docs.css', type: 'text/css'});
		 addTag('link', {rel: 'stylesheet', href: 'css/animations.css', type: 'text/css'});
		 if (jQuery) addTag('script', {src: (debug ? 'components/jquery.js' : 'components/jquery.min.js')});
		 addTag('script', {src: path('angular.js')}, sync);
		 addTag('script', {src: path('angular-resource.js') }, sync);
		 addTag('script', {src: path('angular-route.js') }, sync);

		 */

	function addTag(name, attributes, sync) {
		var el = document.createElement(name),
			attrName;

		for (attrName in attributes) {
			el.setAttribute(attrName, attributes[attrName]);
		}

		sync ? document.write(outerHTML(el)) : headEl.appendChild(el);
	}

	function outerHTML(node){
		// if IE, Chrome take the internal method otherwise build one
		return node.outerHTML || (
			function(n){
				var div = document.createElement('div'), h;
				div.appendChild(n);
				h = div.innerHTML;
				div = null;
				return h;
			})(node);
	}

});