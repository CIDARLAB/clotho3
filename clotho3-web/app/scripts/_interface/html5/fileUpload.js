//note - can't use progress or cancel without extension (shim provided by bower_component)
angular.module('clotho.interface')
.directive('clothoFileUpload', function ($upload, $window, $timeout) {
		return {
		restrict: 'A',
		templateUrl: 'views/_interface/fileUpload.html',
		controller: function ($scope, $element, $attrs) {
			$scope.uploadRightAway = false;

			$scope.upload = [];
			$scope.uploadResult = [];
			$scope.selectedFiles = [];
			$scope.started = [];

			$scope.onFileSelect = function($files) {
				console.log('files changed');
				$scope.selectedFiles = $files;
				$scope.dataUrls = [];
				$scope.inputFiles = [];
				for ( var i = 0; i < $files.length; i++) {
					var $file = $files[i];
					if ($window.FileReader) {

						console.log($file.type);

						var fileReader;
						if ($file.type.indexOf('text') > -1) {
							fileReader = new FileReader();
							function readContents(fileReader, index) {
								fileReader.readAsText($files[i]);
								fileReader.onload = function(e) {
									$timeout(function() {
										$scope.inputFiles[index] = e.target.result;
									});
								}
							}
							readContents(fileReader, i);
						}
						else if ($file.type.indexOf('image') > -1) {
							fileReader = new FileReader();
							fileReader.readAsDataURL($files[i]);
							function setPreview(fileReader, index) {
								fileReader.onload = function(e) {
									$timeout(function() {
										$scope.dataUrls[index] = e.target.result;
									});
								}
							}
							setPreview(fileReader, i);
						}
					}
					$scope.started[i] = false;

					if ($scope.uploadRightAway) {
						$scope.start(i);
					}
				}
			};

			$scope.start = function(index) {
				console.log('starting upload');
				$scope.started[index] = true;
				$scope.upload[index] = $upload.upload({
					url : 'upload/url',
					headers: {
						'myHeaderKey': 'myHeaderVal'
					},
					data : {
						myModel : $scope.fileName
					},
					/*
					//see https://github.com/danialfarid/angular-file-upload/pull/40
					formDataAppender: function(fd, key, val) {
						if (angular.isArray(val)) {
					    angular.forEach(val, function(v) {
					      fd.append(key, v);
					    });
						} else {
							fd.append(key, val);
						}
					},
					*/
					file: $scope.selectedFiles[index],
					fileFormDataName: 'myFile'
				}).then(function(response) {
						$scope.uploadResult.push(response.data.result);
					},
					angular.noop
				)
			}
		},
		link: function (scope, element, attrs) {}
	}
});