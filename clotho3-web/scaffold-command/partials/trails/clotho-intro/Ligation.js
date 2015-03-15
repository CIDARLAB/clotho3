'use strict';

$clotho.extensions.controller('clothoIntro_LigationCtrl', function ($scope) {

	//copied from Digest module for when stripped from client
	$scope.cutMarks = {
		'blunt': {
			mark: '|',
			type: 'Blunt',
			description: 'A cut that results in no "sticky" ends, or overhangs. '
		},
		'main': {
			mark: '^',
			type: 'Main Strand',
			description: 'Denotes a cut on the 5\' -> 3\' (primary) strand, usually visualized as the "top" strand.'
		},
		'comp': {
			mark: '_',
			type: 'Complementary Strand',
			description: 'Denotes a cut on the 3\' -> 5\' (complementary) strand, usually visualized as the "bottom" strand.'
		}
	};
});