angular.module('clotho.interface')
  .directive('darkBackground', function ($document) {

		var body = $document.find('body');

    return {
      restrict: 'A',
      link: function postLink(scope, element, attrs) {
	      body.toggleClass('dark', true);

	      scope.$on("$destroy", function() {
		      body.toggleClass('dark', false);
	      })
      }
    };
  });
