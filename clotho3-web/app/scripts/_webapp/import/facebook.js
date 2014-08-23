'use strict';

angular.module('clotho.webapp')
.controller('FacebookImportCtrl', function ($scope, Clotho, $window) {

		$window.fbAsyncInit = function() {
			FB.init({
				appId      : '1456992994565119',
				xfbml      : true,
				version    : 'v2.0',
				cookie     : true,
				status     : true
			});

			FB.getLoginStatus(function(statusReponse) {
				console.log('statusChangeCallback');
				console.log(statusReponse);
				// The response object is returned with a status field that lets the
				// app know the current login status of the person.
				// Full docs on the response object can be found in the documentation
				// for FB.getLoginStatus().
				if (statusReponse.status === 'connected') {
					// Logged into your app and Facebook.
					console.log('Welcome!  Fetching your information.... ');
					FB.api('/me', function(response) {
						console.log('Successful login: ', response);
						$scope.response = 'welcome ' + response.name;
					});
				} else if (statusReponse.status === 'not_authorized') {
					// The person is logged into Facebook, but not your app.
					$scope.response = 'You are logged into facebook, not this app';
				} else {
					// The person is not logged into Facebook, so we're not sure if
					// they are logged into this app or not.
					$scope.response = 'Not logged into facebook';
				}
			});

			$scope.login = function (cb) {
				FB.login(function (response) {
					console.log('login response', response);
					cb(response);
				}, {scope: 'public_profile email'});
			};

			//log out of facebook entirely
			$scope.logout = FB.logout;
		};



		(function(d, s, id){
			var js, fjs = d.getElementsByTagName(s)[0];
			if (d.getElementById(id)) {return;}
			js = d.createElement(s); js.id = id;
			js.src = "//connect.facebook.net/en_US/sdk.js";
			fjs.parentNode.insertBefore(js, fjs);
		}(document, 'script', 'facebook-jssdk'));

});
