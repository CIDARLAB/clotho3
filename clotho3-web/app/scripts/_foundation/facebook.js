angular.module('clotho.foundation')
.service('Facebook', function ($q, $window, $filter) {

	//todo - incorporate Auth service

	var self = this;

	//future - make provider and allow this to be configurable
	var appId = 1456992994565119;

	var loggingIn = false, 
			loggedIn = false;
		
	var sdkReadyDeferred = $q.defer(),
			sdkReady = sdkReadyDeferred.promise;

	var user,
			userReadyDeferred,
			userReady;

	function setupEmptyUser () {
		user = {};
		userReadyDeferred = $q.defer();
		userReady = userReadyDeferred.promise;
	}

	//init
	setupEmptyUser();

	function fetchMe () {
		loggingIn = true;
		var deferred = $q.defer();

		sdkReady.then(function () {
			FB.api('/me', function (response) {
				handleRetrieveUser(response);
				deferred.resolve(response);
			});
		});

		return deferred.promise;
	}

	//called on FB.login() or if already logged in
	function handleRetrieveUser (response) {
		loggingIn = true;
		loggedIn = false;

		angular.extend(user, response, {
			provider : 'facebook',
			icon : '//graph.facebook.com/' + response.id + '/picture?width=500&height=500'
		});

		userReadyDeferred.resolve(user);
	}

	function handleLoginError (response) {
		loggingIn = false;
		loggedIn = false;

		userReadyDeferred.reject(response);
	}

	//todo - handle logout in session

	//callback for JS API init
	$window.fbAsyncInit = function() {
		FB.init({
			appId      : appId,
			xfbml      : true,
			version    : 'v2.0',
			cookie     : true,
			status     : true
		});

		FB.getLoginStatus(function(statusReponse) {
			//if we want the session id, token or anything... thats in here
			//console.log(statusReponse.authResponse);

			if (statusReponse.status === 'connected') {
				fetchMe();
			} else if (statusReponse.status === 'not_authorized') {
				// The person is logged into Facebook, but not your app.
			} else {
				// The person is not logged into Facebook
			}

			//put this in the callback so we don't call until doing the first check
			sdkReadyDeferred.resolve(FB);
		});

		FB.Event.subscribe('auth.login', function (evt) {
			//todo
		});
		FB.Event.subscribe('auth.logout', function (evt) {
			//todo
		});
};

	//note - calling this without a button click will probably hide the popup
	self.login = function (options) {

		var loginOpts = angular.extend({
			scope: 'public_profile email'
		}, options);

		var loginDeferred = $q.defer();

		sdkReady.then(function (FB) {
			if (loggingIn) {
				userReady.then(function (user) {
					loginDeferred.resolve(user);
				});
			} else {
				FB.login(function (response) {
					console.log('FB login response', response);
					if (response.authResponse) {
						fetchMe().then(function () {
							userReady.then(function(userInfo) {
								loginDeferred.resolve(userInfo);
							});
						});
					} else {
						loginDeferred.reject(response);
					}
				}, loginOpts);
			}
		});

		return loginDeferred.promise;
	};

	self.logout = function (cb) {
		return sdkReady.then(function (FB) {

			FB.getLoginStatus(function (statusReponse) {
				if (statusReponse.status === 'connected') {
					FB.logout(cb)
				}
			});

			//note - likely need to reload the page to clear the cookie etc.

			setupEmptyUser();
		});
	};

	/* The big daddy. Call this to get things started:

	Facebook.getUser()
	.then(function (user) {
		//they are logged in and have authorized app before
	}, function (err) {
		//you will want to log them in
	});

	It is recommended you use both callbacks

	*/
	self.getUser = function () {
		var deferred = $q.defer();
		if (angular.isEmpty(user)) {
			sdkReady.then(function () {
				if (loggingIn) {
					userReady.then(function (user) {
						deferred.resolve(user);
					});
				} else {
					deferred.reject();
				}
			});
		} else {
			deferred.resolve(user);
		}
		return deferred.promise;
	};


	self.convertToPersonSharable = function (input) {
		input = input || user;
		if (angular.isEmpty(user)) {
			return;
		}

		return {
			schema: 'org.clothocad.model.Person',
			id : input.email,
			name : input.name,
			emailAddress : input.email,
			dateCreated : $filter('date')(Date.now().valueOf(), 'yyyy-MM-dd'),
			icon : input.icon,
			social : {
				facebook : input.link
			}
		};
	};

	//IIFE after callbacks etc. set up, will only run when service is instantiated
	(function(d, s, id){
		var js, fjs = d.getElementsByTagName(s)[0];
		if (d.getElementById(id)) {return;}
		js = d.createElement(s); js.id = id;
		js.src = "//connect.facebook.net/en_US/sdk.js";
		fjs.parentNode.insertBefore(js, fjs);
	}(document, 'script', 'facebook-jssdk'));

});