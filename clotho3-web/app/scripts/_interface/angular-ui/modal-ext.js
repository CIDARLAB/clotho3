angular.module('ui.bootstrap-decorate', ['ui.bootstrap'])
	.config(function($provide) {
	$provide.decorator("$modal", ['$delegate', function($delegate) {

		$delegate.messageBox = function(title, message, buttons){
			return $delegate.open({
				templateUrl: 'views/_interface/ui-custom/dialogMessagebox.html',
				controller: 'MessageBoxController',
				resolve: {
					model: function() {
						return {
							title: title,
							message: message,
							buttons: buttons
						};
					}
				}
			});
		};

		$delegate.login = function() {
			return $delegate.open({
				backdrop: true,
				backdropFade: true,
				keyboard: true,
				backdropClick: true,
				templateUrl:  'views/_interface/ui-custom/dialogLogin.html',
				controller: 'DialogLoginController'
			});
		};

		$delegate.serverAlert = function(message) {
			return $delegate.open({
				backdrop: true,
				backdropFade: true,
				keyboard: true,
				backdropClick: true,
				templateUrl: 'views/_interface/ui-custom/dialogMessagebox.html',
				controller: 'ServerAlertController',
				resolve: {
					model: function() {
						return {
							title: "Server Message",
							message: message,
							buttons: [{result:'ok', label: 'OK', cssClass: 'btn-primary'}]
						};
					}
				}
			});
		};

		$delegate.share = function(url) {
			return $delegate.open({
				backdrop: true,
				backdropFade: true,
				keyboard: true,
				backdropClick: true,
				templateUrl:  'views/_interface/ui-custom/dialogShare.html',
				controller: 'DialogShareController',
				resolve: {
					model: function() {
						return {
							url : url
						}
					}
				}
			});
		};

		$delegate.video = function(videoId, videoParams) {
			//todo: preserve aspect ratio
			angular.extend(videoParams, {width: "560"});

			return $delegate.open({
				backdrop: true,
				backdropFade: true,
				keyboard: true,
				backdropClick: true,
				template:  '<div youtube="{{ videoId }}" params="videoParams""></div>',
				controller: 'VideoDialogController',
				resolve: {
					model: function() {
						return {
							videoId : videoId,
							videoParams : videoParams
						}
					}
				}
			});
		};

		$delegate.clothoEdit = function (uuid) {
			return $delegate.open({
				backdrop: true,
				backdropFade: true,
				keyboard: true,
				backdropClick: false,
				templateUrl:  'views/_interface/ui-custom/clothoEditModal.html',
				controller: 'DialogEditController',
				resolve: {
					model: function() {
						return {
							uuid : uuid
						}
					}
				}
			});
		};

		return $delegate;
	}])
});


angular.module('clotho.interface').controller('MessageBoxController', function($scope, $modalInstance, model){
	$scope.title = model.title;
	$scope.message = model.message;
	$scope.buttons = model.buttons;
	$scope.close = function(res){
		$modalInstance.close(res);
	};
});

angular.module('clotho.interface').controller('DialogLoginController', function($scope, $modalInstance, Clotho){
	$scope.close = function(res){
		$modalInstance.close(res);
	};

	$scope.notification = {};
	$scope.cred = {username : "", password: ""};

	$scope.login = function() {
		Clotho.login($scope.cred.username, $scope.cred.password).then(function (result) {
			console.log('run login');
			if (!!result) {
				$scope.notification = {class : "alert-success", message: "Log in Success"};
				dialog.close($scope.cred.username);
			} else {
				$scope.notification = {class : "alert-danger", message: "Log in Error"};
				$scope.cred = {username : "", password: ""};
			}
		});
	};
});

angular.module('clotho.interface').controller('ServerAlertController', function($scope, $modalInstance, model, Clotho){
	$scope.title = model.title;
	$scope.message = model.message;
	$scope.buttons = model.buttons;
	$scope.close = function(res){
		$modalInstance.close(res);
	};

	//todo - more intelligent handling??
	Clotho.listen('serverAlert', function() {
		$scope.close('Another alert appeared');
		Clotho.say($scope.message);
	}, $scope);

});

angular.module('clotho.interface').controller('DialogShareController', function($scope, $modalInstance, model, $location){
	$scope.close = function(result){
		$modalInstance.close(result);
	};

	$scope.customUrl = (model.url && model.url != '') ? model.url : false;

	$scope.social = [
		{
			"name" : "facebook",
			"prefix" : "http://www.facebook.com/sharer.php?u="
		},
		{
			"name" : "google",
			"prefix" : "https://plus.google.com/share?url="
		},
		{
			"name" : "twitter",
			"prefix" : "http://twitter.com/share?url="
		},
		{
			"name" : "linkedin",
			"prefix" : "http://www.linkedin.com/shareArticle?mini=true&url="
		},
		{
			"name" : "digg",
			"prefix" : "http://www.digg.com/submit?url="
		},
		{
			"name" : "reddit",
			"prefix" : "http://reddit.com/submit?url="
		},
		{
			"name" : "email",
			"prefix" : "mailto:?Body="
		}
	];

	$scope.share = function (site) {
		var url = $scope.customUrl ? $scope.customUrl : site.prefix + $location.absUrl();

		$scope.close();

		window.open(url, (site.name == 'email' ? '_self' : "_blank") );
	}

});


angular.module('clotho.interface').controller('VideoDialogController', function($scope, $modalInstance, model){
	$scope.videoId = model.videoId;
	$scope.videoParams = model.videoParams;
	$scope.close = function(res){
		$modalInstance.close(res);
	};
});

angular.module('clotho.interface').controller('DialogEditController', function($scope, $modalInstance, model){
	$scope.uuid = model.uuid;
});


/*
todo - try including in $modal $delegate directly and transcluding

If pass an ID, clotho-show of that id. If no id, then just a modal for transluded content
simply pass true to 'open' and will open a modal right off the bat
 */

angular.module('clotho.interface')
	.directive('clothoModal', function ($modal, $parse) {
		return {
			restrict : 'E',
			transclude : 'true',
			scope: {
				id : '@?',
				open : '@',
				onClose : '=?',
				onOpen : '=?',
				title : '@?',
				model : '@?'
			},
			link: function (scope, element, attrs, nullCtrl, transclude) {

				var modal;

				scope.$watch(function () {
						return $parse(attrs.open);
					}, function (newval, oldval) {
						if (newval != oldval && !!newval) {
							modal = $modal.open({
								backdrop: true,
								backdropFade: true,
								keyboard: true,
								backdropClick: true,
								templateUrl: 'views/_interface/ui-custom/dialogClothoModal.html',
								controller: 'clothoModalController',
								resolve: {
									model: function () {
										return {
											showId : scope.id || '',
											title : scope.title || '',
											model : scope.model | {}
										}
									}
								}
							});

							modal.opened
							.then(function (result) {
								angular.isFunction(scope.onOpen) && scope.onOpen(result);
							});

							modal.result
							.then(function (result) {
								angular.isFunction(scope.onClose) && scope.onClose(result);
							});
						}
					});
			}
		}
	})
	.controller('clothoModalController', function($scope, $modalInstance, model, $compile){

		$scope.title = model.title;
		$scope.model = model.model;

		//if show Id present, use a clotho-show
		$scope.showId = model.showId;
	});