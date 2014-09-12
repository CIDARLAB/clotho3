angular.module('clotho.commandbar')
	.controller('loginCtrl', function ($scope, $timeout, $location, Clotho, ClothoAuth, PubSub, Facebook) {

		$scope.createMode = false;
		$scope.notification = {};

		function resetPassword() {
			$scope.cred.password = '';
			$scope.cred.confirm = '';
		}

		function resetCred() {
			$scope.cred = {username: "", password: "", confirm: "", personId : ""};
		}
		resetCred();

		PubSub.on('auth:login', function () {
			$location.replace();
			$location.path('/settings');
		});

		/*
		//todo - might want to move this to run() clause for app so done immediately
		//todo - work with Facebook + Auth to check on load if someone is logged in

		//when controller is instantiated, try to get info from facebook
		Facebook.getUser().then(function (user) {
			//they have logged into our app before
			$scope.cred = {
				username : user.email
			};
			$scope.facebookRetrieved = true;
		}, function (err) {
			//not logged into our app, show login button
			$scope.showFacebookLogin = true;
		});
		*/

		$scope.login = function () {
			//will only get here if form is valid
			Clotho.login($scope.cred.username, $scope.cred.password)
			.then(function (result) {
				if (result) {
					$scope.notification = {class: "alert-success", message: "Log in Success"};
				} else {
					$scope.notification = {class: "alert-danger", message: "Log in Error"};
					resetPassword();
				}
			}, function (err) {
				//error communicating
			});
		};

		$scope.$watch('cred.username', function (newval) {
			newval && Clotho.get(newval, {mute : true}).then(function (retrieved) {
				if (retrieved && retrieved.id) {
					$scope.retrieved = retrieved;
					$scope.cred.personId = retrieved.id;
					//in case they were in this mode
					$scope.createMode = false;
				} else {
					$scope.retrieved = null;
					$scope.cred.personId = '';
				}
			}, function () {
				$scope.retrieved = null;
				$scope.cred.personId = '';
			});
		});

		$scope.logout = function () {
			Facebook.logout().then(function () {
				resetCred();
				$scope.createMode = false;
				$scope.notification = {};
				$scope.showFacebookLogin = true;
				$scope.retrieved = null;
				$scope.facebookRetrieved = false;
			});
		};


		//todo - combine with facebook import
		$scope.facebookLogin = function () {
			Facebook.login().then(function (user) {

				$scope.cred.username = user.email;
				$scope.facebookRetrieved = true;

				//check if they exist in clotho
				Clotho.get(user.email, {mute : true}).then(function (retrieved) {

					if (!angular.isEmpty(retrieved)) {
						console.log('get returned', retrieved);

						$scope.notification = {
							class: "alert-success",
							message: "Enter your password to login!"
						};

						$scope.cred.username = user.email;
						$scope.cred.personId = user.email;
					} else {
						console.log('get was null, creating');

						var person = Facebook.convertToPersonSharable(user);

						Clotho.create(person)
							.then(function (id) {
								$scope.notification = {
									class: "alert-success",
									message: "Account added to Clotho! Choose a password"
								};
								$scope.retrieved = person;
								$scope.cred.username = user.email;
								$scope.cred.personId = id;
							}, function () {
								//already existed, so just add info to form
								$scope.notification = {
									class: "alert-success",
									message: "You've imported that account already! Choose a password"
								};
								$scope.retrieved = person;
								$scope.cred.username = user.email;
								$scope.cred.personId = person.id;
							})
							.then(function () {
								$scope.createMode = true;
							});
					}
				});
			}, function (err) {
				$scope.notification = {
					class: "alert-danger",
					message: "Logging into facebook went wrong..."
				};
				//todo - we will need to re-request if denied: see the API
			});
		};

		$scope.createAccount = function () {
      //should only get here if form is valid

      /*
      //todo - incorporate associated user into flow + errors
      //check associated person exists
      var hasAssociated = !angular.isEmpty($scope.cred.personId);

      if (hasAssociated) {
        Clotho.get($scope.cred.personId)
        .then(function (retrieved) {
          //double check that associated person exists
          hasAssociated = !!retrieved;
        });
      }
      */

      //try to create
      Clotho.createUser($scope.cred.username, $scope.cred.password)
        .then(function (response) {

          console.log('create user?', response);

          //todo

          if (response) {
            $scope.notification = {
              class: "alert-success",
              message: "User " + $scope.cred.username + "created!"
            };
          } else {
            $scope.notification = {
              class: "alert-error",
              message: "Account creation unsuccessful"
            };
          }

        }, function (err) {
          $scope.notification = {
            class: "alert-error",
            message: "Error Creating... check console"
          };
          console.error('account creation error', err);
        });
		};

	});
