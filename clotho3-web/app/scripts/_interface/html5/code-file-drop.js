//note - can't use progress or cancel without extension (shim provided by bower_component)
angular.module('clotho.interface')
	.directive('functionCodeDrop', function ($upload, $window, $timeout) {
		return {
			restrict: 'A',
			templateUrl: '../../../views/_interface/codeDrop.html',
			scope : {
				updateOnRead : '=',
				showDrop: '@',
				showBrowser: '@'
			},
			controller: function ($scope, $element, $attrs) {
				$scope.uploadRightAway = false;

				$scope.selectedFiles = [];
				$scope.inputFiles = [];

				function readContents(fileReader, file) {
					fileReader.readAsText(file);
					fileReader.onload = function(e) {
						$timeout(function() {
							$scope.updateOnRead = e.target.result;
						});
					}
				}

				$scope.onFileSelect = function($files) {
					$scope.selectedFiles = $files;

					for ( var i = 0; i < $files.length; i++) {
						var $file = $files[i];
						if ($window.FileReader && $file.type.indexOf('text') > -1) {
							var fileReader = new FileReader();
							readContents(fileReader, $files[i]);
						}
					}
				};
			},
			link: function (scope, element, attrs) {}
		}
	});