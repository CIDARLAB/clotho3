//note : localStorage only supports strings, so we need to manually serialize and deserialize

angular.module('clotho.core').service('clothoLocalStorage',
	function($window, PubSub, Debug, $document) {

		return ($window.$clotho.$localStorage) ? $window.$clotho.$localStorage : $window.$clotho.$localStorage = generateLocalStorage();

		function generateLocalStorage() {

			var Debugger = new Debug('clothoLocalStorage', '#88cc88');

			//defaults
			var refStorage = $window.localStorage;
			var serializer = JSON;
			//prefix to prevent name-clashes
			var prefix = "clotho_";


			function broadcastModelUpdate(uuid, obj) {
				PubSub.trigger("update", uuid);
				PubSub.trigger("update:" + uuid, angular.copy(obj));
			}

			// --- Check for support ---
			// FUTURE - add this in later to be more robust
			// e.g. https://github.com/grevory/angular-local-storage/blob/master/localStorageModule.js
			// and could fallback to cookies for small strings (<4kb)

			/*
			Checks the browser to see if local storage is supported. Adopted from Modernizer.

			Note that this will trigger cross browser events that we pick up in this service.

			  In FF4, if disabled, window.localStorage should === null.

				Normally, we could not test that directly and need to do a
				 `('localStorage' in window) && ` test first because otherwise Firefox will
				 throw bugzil.la/365772 if cookies are disabled

				Also in iOS5 Private Browsing mode, attempting to use localStorage.setItem
				will throw the exception:
				 QUOTA_EXCEEDED_ERRROR DOM Exception 22.
				Peculiarly, getItem and removeItem calls do not throw.

				Because we are forced to try/catch this, we'll go aggressive.

				Just FWIW: IE8 Compat mode supports these features completely:
				 www.quirksmode.org/dom/html5.html
				But IE8 doesn't support either with local files
			*/
			var browserSupportsLocalStorage = function () {
				try{
					refStorage.setItem('test', '7');
					if(refStorage.getItem('test')=== '7'){
						refStorage.removeItem('test');
						return true;
					}
				}
				catch(er){}
				return false;
			};

			// Checks the browser to see if cookies are supported
			// look into angular $cookies if we want to implement this as a backup
			// note the limitations of cookies before committing to that
			var browserSupportsCookies = function () {
				try {
					return navigator.cookieEnabled ||
						("cookie" in $document && ($document.cookie.length > 0 ||
							($document.cookie = "test").indexOf.call($document.cookie, "test") > -1));
				} catch (e) {
					return false;
				}
			};

			Debugger.log("localStorage support? " + browserSupportsLocalStorage());

			// --- local storage interface ----

			//remove items with matching prefix, return true if item(s) removed
			var clear = function () {
				var flag = false;
				for ( var i = 0, len = refStorage.length; i < len; ++i ) {
					var key = refStorage.key(i);
					if (angular.isString(key) && key.substring(0, prefix.length) == prefix) {
						refStorage.removeItem(key);
						--i; --len;
						flag = true;
					}
				}
				return flag;
			};

			// returns an item, or optional defaultValue if not found
			var getItem = function (key, defaultValue) {
				if (angular.isEmpty(key)) {
					return null;
				}
				var value = refStorage.getItem(prefix + key);
				if (angular.isEmpty(value)) {
					return(angular.isDefined(defaultValue) ? defaultValue : null );
				} else {
					return angular.isObject(value) ? value : serializer.parse(value);
				}
			};

			// return a boolean whether a given object exists
			var hasItem = function (key) {
				if (angular.isEmpty(key)) {
					return false;
				}
				var temp = refStorage.getItem(prefix + key);
				return ( angular.isDefined(temp) && !angular.isEmpty(temp) );
			};

			//remove a given item from storage
			var removeItem = function (key) {
				if (angular.isEmpty(key)) {
					return false;
				}
				refStorage.removeItem(prefix + key);
				return true;
			};

			//adds an item to storage, automatically serializing.
			//NOTE -  cannot add functions or private variables
			var setItem = function (key, value) {
				//prevent adding of empty objects
				if (angular.isEmpty(key)) {
					return false;
				}
				if (angular.isEmpty(value)) {
					return false;
				}
				value = angular.isString(value) ? value : serializer.stringify(value);
				try {
					refStorage.setItem(prefix + key, value);
				} catch (err) {
					Debugger.error('couldnt save ' + key, err);
				}
				return true;
			};

			// ----- LOCAL STORAGE EVENTS ----
			/* NOTES
			 this half is for model update broadcasting across pages (via localStorage)
			 the other half is using PUB/SUB (within pages + to bubble up updates)

			 note - implementation of localStorage events varies & is unreliable across browsers
			 - event handlers only invoked for current window / tab where data is written / deleted (even though spec states all windows)
			 - may not be called when a key is updated, and only when it is created / deleted
			 - storage deletion only returns key, not deleted value
			 */

			var handle_storage_change = function (e) {
				if (!e) {
					e = $window.event;
				}

				// avoid non-prefixed items
				if (e.key.indexOf(prefix) < 0) {
					return;
				}

				var uuid = e.key.replace(prefix, '') || '';
				//note - this relies on the object existing in localStorage (which it won't if reset on page load)
				var obj = getItem(uuid);
				Debugger.log("handle_storage_event for " + uuid, obj);
				broadcastModelUpdate(uuid, obj)
			};

			if ($window.addEventListener) {
				$window.addEventListener("storage", handle_storage_change, false);
			} else {
				$window.attachEvent("onstorage", handle_storage_change); //IE8
			}


			return {
				getPrefix: function () {
					return prefix
				},
				isSupported: browserSupportsLocalStorage,
				clear: clear,
				getItem: getItem,
				hasItem: hasItem,
				removeItem: removeItem,
				setItem: setItem
			}
		}
});
