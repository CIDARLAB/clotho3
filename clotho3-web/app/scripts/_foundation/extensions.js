/*!
 $script.js v1.3 -- Async script loader
 usage: https://github.com/ded/script.js
 */
!function(a,b,c){function t(a,c){var e=b.createElement("script"),f=j;e.onload=e.onerror=e[o]=function(){e[m]&&!/^c|loade/.test(e[m])||f||(e.onload=e[o]=null,f=1,c())},e.async=1,e.src=a,d.insertBefore(e,d.firstChild)}function q(a,b){p(a,function(a){return!b(a)})}var d=b.getElementsByTagName("head")[0],e={},f={},g={},h={},i="string",j=!1,k="push",l="DOMContentLoaded",m="readyState",n="addEventListener",o="onreadystatechange",p=function(a,b){for(var c=0,d=a.length;c<d;++c)if(!b(a[c]))return j;return 1};!b[m]&&b[n]&&(b[n](l,function r(){b.removeEventListener(l,r,j),b[m]="complete"},j),b[m]="loading");var s=function(a,b,d){function o(){if(!--m){e[l]=1,j&&j();for(var a in g)p(a.split("|"),n)&&!q(g[a],n)&&(g[a]=[])}}function n(a){return a.call?a():e[a]}a=a[k]?a:[a];var i=b&&b.call,j=i?b:d,l=i?a.join(""):b,m=a.length;c(function(){q(a,function(a){h[a]?(l&&(f[l]=1),o()):(h[a]=1,l&&(f[l]=1),t(s.path?s.path+a+".js":a,o))})},0);return s};s.get=t,s.ready=function(a,b,c){a=a[k]?a:[a];var d=[];!q(a,function(a){e[a]||d[k](a)})&&p(a,function(a){return e[a]})?b():!function(a){g[a]=g[a]||[],g[a][k](b),c&&c(d)}(a.join("|"));return s};var u=a.$script;s.noConflict=function(){a.$script=u;return this},typeof module!="undefined"&&module.exports?module.exports=s:a.$script=s}(this,document,setTimeout);

/**
 * @ngdoc module
 * @name clotho.extensions
 * @description
 * Module for lazy-registering angular components to a web-app.
 *
 * note - this module should not contain any components (controller, directives, etc.) to be loaded at initial bootstrap. $inject calls are re-routed through provider to add them manually
 *
 * note - you can always add a new module when you bootstrap a new widget. That widget can inherit the functionality of the clotho webapp by including the appropriate package (at time of writing, clotho.fullPackage). Then you can have config and run blocks, ignore routing, run autonomously, configure providers, etc.
 *
 */
angular.module('clotho.extensions', [])
.config(function($controllerProvider, $compileProvider, $filterProvider, $provide, $injector, $animateProvider) {

	window.$clotho.extensions = $clotho.extensions = {};

	$clotho.extensions.providers = {
		$controllerProvider: $controllerProvider,
		$compileProvider: $compileProvider,
		$filterProvider: $filterProvider,
		$provide: $provide,
    $injector: $injector,
    $animateProvider: $animateProvider
	};

	// adapted from Ifeanyi Isitor: http://ify.io/lazy-loading-in-angularjs/
	// preserve the original references. angular.module('clotho.extensions').controller() is the same as $clotho.extensions._controller(). You can use these references and process the queue manually. Recommended you use the provider-based interface below though, and download components which take the form $clotho.extensions.controller() [without the underscore] which will be added automatically
	$clotho.extensions._controller = $clotho.extensions.controller;
	$clotho.extensions._service = $clotho.extensions.service;
	$clotho.extensions._factory = $clotho.extensions.factory;
	$clotho.extensions._value = $clotho.extensions.value;
	$clotho.extensions._directive = $clotho.extensions.directive;
	$clotho.extensions._animation = $clotho.extensions.animation;

	// Convert controllers, services, factories, directives, filters, animations to provider-based alternatives, so that when they are loaded later, they are automatically registered
	$clotho.extensions.controller = $controllerProvider.register;
	$clotho.extensions.service = $provide.service;
	$clotho.extensions.factory = $provide.factory;
	$clotho.extensions.value = $provide.value;
	$clotho.extensions.directive = $compileProvider.directive;
	$clotho.extensions.filter = $filterProvider.filter;
	$clotho.extensions.animation = $filterProvider.animation;

})
.run(function($rootScope, $q, $timeout, $templateCache, $http, $cacheFactory, $rootElement) {

	var registeredQueueLength = 0,
      timeoutTime = 5000,
      fileCache = $cacheFactory('fileUrls');

  //todo - add service

  //todo - mainModule._configBlocks for angular 1.3 (see ocLazyLoad)

	//these functions are for mixing in components that are not added via the provider interface above
  function getQueue () {
		return angular.module('clotho.extensions')._invokeQueue;
	}

  function updateQueueLength () {
    registeredQueueLength = getQueue().length;
  }

  //init
  updateQueueLength();

	//Only needs to be called when not using the provider interface above, i.e. components that have been prefaced, e.g. angular.module('clotho.extensions').controller() instead of $clotho.extensions.controller()
  //todo - verify that processing queue after adding components using provider-based mechanism doesn't cause problems with queue length mismatch
	var processQueue = function() {
		var queue = getQueue();

		for(var i=registeredQueueLength;i<queue.length;i++) {
			var call = queue[i];
			// call is in the form [providerName, providerFunc, providerArguments]
			var provider = $clotho.extensions.providers[call[0]];
			if(provider) {
				// e.g. $controllerProvider.register("Ctrl", function() { ... })
				//console.log("provider exists");
				provider[call[1]].apply(provider, call[2]);
			}
		}

    updateQueueLength();
	};

  function cacheBust(url) {
    var dc = new Date().getTime();
    if(url.indexOf('?') >= 0) {
      if(url.substring(0, url.length - 1) === '&') {
        return url + '_dc=' + dc;
      }
      return url + '&_dc=' + dc;
    } else {
      return url + '?_dc=' + dc;
    }
  }

  function detFileType (url) {
    //The extension is always the last characters before the ? and after a period.
    //accounting for the possibility of a period in the query string
    var b = url.split('?')[0];
    return b.substr(b.lastIndexOf('.')+1);
  }

  function mixinJs (url, callback) {
    $script(url, function () {
      callback(url);
    });
  }

  //unlike mixin, scripts are timestamped so will be run every time
  function loadJs (url, callback) {
    mixinJs(cacheBust(url), callback);
  }

  function loadCss (url, callback) {
    if (document.createStyleSheet) { //IE
      document.createStyleSheet(url);
    } else {
      var link = document.createElement("link");
      link.type = "text/css";
      link.rel = "stylesheet";
      link.href = url;
      document.getElementsByTagName("head")[0].appendChild(link);
    }
    callback(url);
  }

  function loadHtml (url, forceName, callback) {
    $http.get(url)
      .success(function(data) {
        $templateCache.put(forceName || url, data);
        callback(url);
      })
      .error(function(data) {callback(null)});
  }

    /**
     *
     * @param url
     * @param params {object} Supported params
     *
     * noCache (*) {Boolean} -- do not cache resource
     * cacheName (*) {string} -- name to cache under (for de-duping resources)
     * mixin (scripts) {Boolean}
     * templateName (html) {string}
     *
     * @returns {*}
     */
  //handles a single file. returns a promise
  function loadFilePath (url, params) {
    if (angular.isUndefined(url) || url.length == 0) {
      return $q.when(null);
    }

    params = angular.extend({
      mixin: true
    }, params);

    //check cache
    if (params.mixin !== false && fileCache.get(url)) {
      return $q.when(url);
    }

    var deferred = $q.defer(),
        timeoutPromise = $timeout(function() { deferred.reject(null) }, timeoutTime);

    function handleLoadCallback () {
      $timeout.cancel(timeoutPromise);
      (params.noCache !== false) && fileCache.put(params.cacheName || url);
      $rootScope.$evalAsync(deferred.resolve(url)); //future - remove $evalAsync
    }

    switch (detFileType(url)) {
      case 'js' : {
        if (!!params.mixin) {
          mixinJs(url, handleLoadCallback);
        } else {
          loadJs(url, handleLoadCallback);
        }
        break;
      }
      case 'css' : {
        loadCss(url, handleLoadCallback);
        break;
      }
      case 'html' : {
        loadHtml(url, params.templateName, handleLoadCallback);
        break;
      }
      default : {
        deferred.reject(new Error('Dont know how to handle ' + url));
      }
    }

    return deferred.promise;
  }

  //handle multiple files. returns array of promises for each url
  function loadFiles (urls, params) {
    urls = angular.isArray(urls) ? urls : [urls];
    var promises = [];
    angular.forEach(urls, function (url) {
      promises.push(loadFilePath(url, params));
    });
    return $q.all(promises);
  }

	/**
	 * @name downloadDependencies
	 *
	 * @param {Object} dependencies [see format below]
	 *
	 * @returns {Promise} promise which resolves when css, mixins, scripts are downloaded, and a $timeout() for the onload script
	 *
	 * @description
	 *
	 * Given a dependencies object, e.g.: {
	 *   css : <dep>,
	 *   mixin : <dep>,
	 *   script : <dep>,
	 *   onload : <dep>,
	 * }
	 *
	 * where <dep> is either a file path, or array of paths, this function downloads all dependencies, and returns a promise containing a $timeout() which will run the onload() on the next $digest()
	 */
	$clotho.extensions.downloadDependencies = function downloadDependencies(dependencies) {
		if (angular.isObject(dependencies) && Object.keys(dependencies).length > 0) {
      return $q.all([
        $clotho.extensions.css(dependencies.css),
        $clotho.extensions.mixin(dependencies.mixin),
        $clotho.extensions.script(dependencies.script)
      ]).then(function () {
        return function () {
          $timeout(function () {
            $clotho.extensions.script(dependencies.onload);
          });
        };
      });
		} else {
      return loadFiles(dependencies);
    }
	};

	/**
	 * @name $clotho.extensions.mixin
	 *
	 * @description Will download URLs only once
	 *
	 * @param urls URLs of dependencies. Only downloaded if hasn't been already
	 * @returns {Promise} Promise to be fulfilled on successful addition, value is urls passed
	 */
	$clotho.extensions.mixin = loadFiles;

	/**
	 * @name $clotho.extensions.script
	 *
	 * @description Downloads and executes a script or scripts, using cache-busting. Timestamp is appended to the script, to ensure it will run each time.
	 *
	 * @param {string|Array} urls URLs of scripts to be executed
	 * @returns {Promise} Promise which is resolved when all scripts have been executed
	 */
	$clotho.extensions.script = function (urls) {
    return loadFiles(urls, {mixin : false})
  };

	/**
	 * @name $clotho.extensions.css
	 *
	 * @description Downloads and appends a CSS file to the page head, so that it will be applied properly
   *
   * Style sheets will be applied globally.
   * You should use shadowDOM and import manually for style encapsulation.
   *
   * Prefix all styles so they do not affect other components.
	 *
	 * @param {string|Array} urls URL of CSS file
	 *
	 * @returns {Promise} Promise to be fulfilled once CSS files are downloaded and appended
	 */
	$clotho.extensions.css = loadFiles;

	/**
	 * @name $clotho.extensions.cache
	 *
	 * @description Downloads caches an angular template for later use. Forces addition to the cache, under ID of the passed URL. So, to use later, e.g. use <div ng-include="url_you_passed"></div>. Note that the template is cached in the primary app, so to access it in a separately bootstrapped app, you'll need to list the appropriate angular module as a dependency
	 *
	 * @param {URL} url URL of angular template file
   * @param {String} forceName
   * @param {Function} callback
	 *
	 * @returns {Promise} Promise to be fulfilled once CSS files are downloaded and appended
	 */
	$clotho.extensions.cache = function (url, templateName) {
    return loadFiles(url, {templateName : templateName});
  };

  /**
   * @name determineUrlExtension
   * @description Given an arbitrary URL, determine file name extension
   * @param url {String}
   * @returns {String}
   */
	$clotho.extensions.determineUrlExtension = detFileType;

  // (re)compile an element, with optional locals for the scope
  // note - create a new isolate scope (child of this module's $rootElement's $rootScope), deleting scope and children scopes if present
  // note - can't really scrub DOM of compiled code, so should use bootstrap for something that needs to be compiled fresh -- do not double-compile
  $clotho.extensions.recompile = function(element, args) {
    if (angular.isUndefined(element)) {return;}
    args = args || {};

    if (element.scope()) {element.scope().$destroy();}

    $rootElement.injector().invoke(['$compile', '$rootScope', function($compile, $rootScope) {
      var scope = $rootScope.$new(true);
      angular.extend(scope, args);
      $compile(element)(scope);
      $rootScope.$apply();
    }]);
  };

  // If you just want to add a property to the $rootScope, given an object to extend() it
  // not really recommended as will slow the $digest cycle
  $clotho.extensions.extendRootscope = function(args) {
    angular.extend($rootScope, args);
  };
});
