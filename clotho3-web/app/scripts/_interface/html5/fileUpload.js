/**
 * WORK IN PROGRESS
 *
 * @name clothoFileUpload
 * @ngdoc directive
 * @module clotho.interface
 *
 * @dependency ng-file-drop and $upload
 * relies on ng-file-drop directive (bower components)
 *
 * @description
 * Directive to wrap a file drop. Will read file contents, and can upload or pass back to function.
 *
 * To be meaningful, must either pass startFunction OR uploadUrl and uploadHeaders
 *
 * @attr auto-upload-visible
 * if present, show check to upload all selected files immediately
 * @attr replace-on-select
 * if present, will replace pending list when select new files (default is add to list)
 *
 * @attr on-select-callback {Function}
 * Run after selected and before processing. Passed 'files' which were just selected
 * @attr on-select-each-callback {Function}
 * Run for each file after fileReader has processed (so may not be in order). passed 'file' and 'index'
 *
 * @attr start-function {Function}
 * Function to run instead of upload. If provided, file will not be uploaded. Passed 'file' and 'index' and 'content'
 * @attr on-start {Function}
 * Run when start function/upload
 *
 * @attr upload-url {String}
 * URL to upload file to
 * @attr upload-headers {Object}
 * Headers object to use for upload
 * @attr on-upload-success {Function} Function to run on successful upload. passed 'index' 'file' and 'result'
 * @attr on-upload-error {Function} Function to run on errored upload. passed 'index' 'file' and 'result'
 *
 * todo - ability to remove file from list
 *
 * //note - can't use progress or cancel without extension (shim provided by bower_component)
 *
 * //future (perf) - should build in timeouts for large files
 */
angular.module('clotho.interface')
.directive('clothoFileUpload', function ($upload, $window, $timeout) {
		return {
		restrict: 'E',
		templateUrl: 'views/_interface/fileUpload.html',
		scope : {
			onSelectCallback : '&?',
			onSelectEachCallback : '&?',
			onStart : '&?',
			startFunction : '&?',
			uploadUrl : '=?',
			uploadHeaders : '=?',
			onUploadSuccess : '&?',
			onUploadError : '&?'
		},
		controller: function ($scope, $element, $attrs) {
			$scope.uploadRightAway = false;

			function createLists () {
				$scope.upload = [];
				$scope.uploadResult = [];
				$scope.selectedFiles = [];
				$scope.started = [];

				//track whether an index is a dataUrl
				$scope.dataUrls = [];
				//file contents (string for text, base64 for image, binary otherwise)
				$scope.inputFiles = [];
			}

			$scope.onFileSelect = function($files) {

				if (angular.isDefined($scope.onSelectCallback)) {
					$scope.onSelectEachCallback({
						files: $files
					});
				}

				//index for for loop to go through new files
				var i;
				//future - ability to concat files... need to update all lists
				if (angular.isUndefined($scope.selectedFiles) || angular.isDefined($attrs.replaceOnSelect)) {
					createLists();
					$scope.selectedFiles = $files;
					i = 0;
				} else {
					i = $scope.selectedFiles.length;
					$scope.selectedFiles = $scope.selectedFiles.concat($files);
				}

				for (var f = 0; f < $files.length; f++, i++) {
					var $file = $scope.selectedFiles[i];
					if ($window.FileReader) {

						var fileReader = new FileReader();

						function readFileContents (fileReader, file, index) {
							fileReader.onload = function(e) {
								$timeout(function() {
									$scope.inputFiles[index] = e.target.result;
								});
							};

							if (angular.isDefined($scope.onSelectEachCallback)) {
								$scope.onSelectEachCallback({
									index: index,
									file: file
								});
							}
						}

						if ($file.type.indexOf('text') > -1) {
							fileReader.readAsText($file);
							readFileContents(fileReader, $file, i);
						}
						else if ($file.type.indexOf('image') > -1) {
							$scope.dataUrls[i] = true;
							fileReader.readAsDataURL($file);
							readFileContents(fileReader, $file, i);
						}
						else {
							//try to read as text.. handle genbanks
							fileReader.readAsBinaryString($file);
							readFileContents(fileReader, $file, i);
						}
					}
					$scope.started[i] = false;

					if ($scope.uploadRightAway) {
						$scope.start(i);
					}
				}
			};

			function uploadFile (index) {

				if (!angular.isDefined($scope.uploadUrl)) {
					console.warn('file could not be uploaded, no url provided');
					return;
				}

				$scope.upload[index] = $upload.upload({
					url : $scope.uploadUrl,
					headers: $scope.uploadHeaders,
					data : {
						myModel : $scope.selectedFiles[index].name
					},
					file: $scope.selectedFiles[index],
					fileFormDataName: 'fileUpload' + Date.now().valueOf()
				})
				.then(function(response) {
					$scope.uploadResult.push(response.data.result);
					if (angular.isDefined($scope.onUploadSuccess)) {
						$scope.onUploadSuccess({
							index : index,
							result : response.data.result,
							file : $scope.selectedFiles[index]
						});
					}
				}, function (err) {
					console.warn('there was an error uploading', index, $scope.selectedFiles[index]);
					if (angular.isDefined($scope.onUploadError)) {
						$scope.onUploadError({
							index : index,
							result : err,
							file : $scope.selectedFiles[index]
						});
					}
				});
			}

			$scope.start = function(index) {

				//check this hasn't been started already
				if ($scope.started[index]) {
					return;
				}

				$scope.started[index] = true;
				if (angular.isDefined($scope.startFunction)) {
					$scope.startFunction({
						index : index,
						file : $scope.selectedFiles[index],
						content : $scope.inputFiles[index]
					})
				} else {
					uploadFile(index);
				}

				if (angular.isDefined($scope.onStart)) {
					$scope.onStart({
						index : index,
						file : $scope.selectedFiles[index],
						content : $scope.inputFiles[index]
					});
				}
			};

			$scope.processAll = function () {
				angular.forEach($scope.selectedFiles, function (file, index) {
					$scope.start(index);
				});
			};

		},
		link: function (scope, element, attrs) {
			scope.autoUploadVisible = angular.isDefined(attrs.autoUploadVisible);

			scope.previewFile = function(index) {
				scope.showPreview = true;
				scope.modalFile = scope.inputFiles[index];
				scope.modalFileIsImage = !!scope.dataUrls[index];
			};
		}
	}
});