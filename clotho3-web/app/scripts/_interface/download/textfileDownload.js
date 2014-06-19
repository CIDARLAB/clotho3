angular.module('clotho.interface')
/**
 * @ngdoc directive
 * @name textfileDownload
 *
 * @description
 * Given a string as content, and an optional name, creates a download link which will download the content as a file txt file.
 */
.directive('textfileDownload', function ($window, $document) {

	const MIME_TYPE = 'text/plain';

	return {
		restrict : 'E',
		scope : {
			content : '=',
			name : '@?'
		},
		link : function (scope, element, attrs) {

			//simple check for necessary APIs
			if (angular.isUndefined($window.Blob) || angular.isUndefined($window.URL.createObjectURL)) {
				return;
			}

			function createLink () {
				var contentText = scope.content;
				var name = scope.name;
				if (angular.isUndefined(contentText)) {
					return;
				}
				var blobby = new $window.Blob([contentText], {type : MIME_TYPE});
				var link = $document[0].createElement('a');
				link.download = name;
				link.alt = 'Download as File';
				link.href = $window.URL.createObjectURL(blobby);
				link.innerHTML = '<span class="glyphicon glyphicon-download"></span>';
				element.replaceWith(link);
			}

			scope.$watch('content', createLink);
			scope.$watch('name', createLink);
		}
	}
});
