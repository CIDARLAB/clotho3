'use strict';

console.log('loaded');

$clotho.extensions.controller('clothoTool_digestCuts', function($scope, $focus, $timeout, $modal, Digest) {
	$scope.Digest = Digest;

	$scope.demo = {};
	$scope.demo.digestEnz = Digest.enzymes.EcoRI;
	$scope.demo.digestSeq = 'acaacgtctcacggatccagtcggaattctacatgcatcgatcgacggatccagatcgactagc';
});