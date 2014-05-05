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


angular.module('ui.bootstrap-decorate').controller('MessageBoxController', function($scope, $modalInstance, model){
	$scope.title = model.title;
	$scope.message = model.message;
	$scope.buttons = model.buttons;
	$scope.close = function(res){
		$modalInstance.close(res);
	};
});

angular.module('ui.bootstrap-decorate').controller('DialogShareController', function($scope, $modalInstance, model, $location, $window){
	$scope.close = function(result){
		$modalInstance.close(result);
	};
});


angular.module('ui.bootstrap-decorate').controller('VideoDialogController', function($scope, $modalInstance, model){
	$scope.videoId = model.videoId;
	$scope.videoParams = model.videoParams;
	$scope.close = function(res){
		$modalInstance.close(res);
	};
});

angular.module('ui.bootstrap-decorate').controller('DialogEditController', function($scope, $modalInstance, model){
	$scope.uuid = model.uuid;
});