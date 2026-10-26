;;;; SPDX-FileCopyrightText: Copyright (c) 2026 Genworks International
;;;; SPDX-License-Identifier: Apache-2.0
;;;;
;;;; Every NLib routine the Gendl/GDL stack reaches, with its Lisp name, its
;;;; source, its porting wave and status (:done, :core, :vanish, :port).
;;;; Wave 0 :port is the remaining infrastructure; waves 1-6 follow the
;;;; dependency order between NLib files.  Generated; edit :status as work lands.

(
 (:c "N_BasisGetLongestAndShortestSpans" :lisp "basis-get-longest-and-shortest-spans" :file "NL_BasisAdv.cpp" :line 47 :lines 22 :wave 0 :status :done
   :callees ("N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots"))
 (:c "N_GetCompatibleKnotArrayMult" :lisp "get-compatible-knot-array-mult" :file "NL_BasisAdv.cpp" :line 119 :lines 69 :wave 0 :status :done
   :callees ("N_BasisReparam" "N_CreateInterval" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_BasisIncreaseKnotMult" :lisp "basis-increase-knot-mult" :file "NL_BasisAdv.cpp" :line 585 :lines 28 :wave 0 :status :done
   :callees ("N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_BasisEvalArray" :lisp "basis-eval-array" :file "NL_BasisAdv.cpp" :line 964 :lines 75 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_BasisFindSpan" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots"))
 (:c "N_BasisDerivsArray" :lisp "basis-derivs-array" :file "NL_BasisAdv.cpp" :line 1122 :lines 118 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisFindSpan" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SwapIntegers"))
 (:c "N_KnotVectorIsParamOutOfBounds" :lisp "knot-vector-is-param-out-of-bounds" :file "NL_BasisBasic.cpp" :line 55 :lines 13 :wave 0 :status :done
   :callees ("N_ErrSet" "N_KnotVectorGetKnots"))
 (:c "N_KnotVectorIsEndParam" :lisp "knot-vector-is-end-param" :file "NL_BasisBasic.cpp" :line 110 :lines 16 :wave 0 :status :done
   :callees ("N_ErrSet" "N_KnotVectorGetKnots"))
 (:c "N_KnotVectorIsValid" :lisp "knot-vector-is-valid" :file "NL_BasisBasic.cpp" :line 172 :lines 43 :wave 0 :status :done
   :callees ("N_ErrSet" "N_KnotVectorGetKnots"))
 (:c "N_AllocKnotVector" :lisp "alloc-knot-vector" :file "NL_BasisBasic.cpp" :line 260 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocKnotVectorAndArray" :lisp "alloc-knot-vector-and-array" :file "NL_BasisBasic.cpp" :line 324 :lines 13 :wave 0 :status :vanish
   :callees ("N_AllocKnotVector" "N_AllocReal1dArray" "N_KnotVectorFromRealArray"))
 (:c "N_Alloc1dArrayKnotVectors" :lisp "alloc-1d-array-knot-vectors" :file "NL_BasisBasic.cpp" :line 380 :lines 35 :wave 0 :status :vanish
   :callees ("N_AllocKnotVectorAndArray" "N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_Alloc1dArrayKnotVectPtrs" :lisp "alloc-1d-array-knot-vect-ptrs" :file "NL_BasisBasic.cpp" :line 457 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_KnotVectorGetKnots" :lisp "knot-vector-get-knots" :file "NL_BasisBasic.cpp" :line 519 :lines 5 :wave 0 :status :vanish
   :callees ())
 (:c "N_KnotVectorFromRealArray" :lisp "knot-vector-from-real-array" :file "NL_BasisBasic.cpp" :line 557 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_SetKnotIndex" :lisp "set-knot-index" :file "NL_BasisBasic.cpp" :line 592 :lines 4 :wave 0 :status :vanish
   :callees ())
 (:c "N_FreeKnotVector" :lisp "free-knot-vector" :file "NL_BasisBasic.cpp" :line 626 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_KnotsCopy" :lisp "knots-copy" :file "NL_BasisBasic.cpp" :line 763 :lines 15 :wave 0 :status :done
   :callees ("N_KnotVectorGetKnots" "N_KnotsCheck"))
 (:c "N_KnotsCheck" :lisp "knots-check" :file "NL_BasisBasic.cpp" :line 822 :lines 27 :wave 0 :status :done
   :callees ("N_AllocKnotVectorAndArray" "N_ErrSet" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_BasisFindSpan" :lisp "basis-find-span" :file "NL_BasisBasic.cpp" :line 1429 :lines 70 :wave 0 :status :done
   :callees ())
 (:c "N_BasisFindSpanAndMult" :lisp "basis-find-span-and-mult" :file "NL_BasisBasic.cpp" :line 1597 :lines 46 :wave 0 :status :done
   :callees ("N_BasisFindSpan" "N_KnotVectorGetKnots"))
 (:c "N_BasisGetSpanCount" :lisp "basis-get-span-count" :file "NL_BasisBasic.cpp" :line 1696 :lines 12 :wave 0 :status :done
   :callees ("N_KnotVectorGetKnots"))
 (:c "N_BasisFindGlobalMax" :lisp "basis-find-global-max" :file "NL_BasisBasic.cpp" :line 1751 :lines 94 :wave 0 :status :done
   :callees ("N_BasisIDerivs" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_KnotVectorGetKnots"))
 (:c "N_BasisFindAllSpanMaxima" :lisp "basis-find-all-span-maxima" :file "NL_BasisBasic.cpp" :line 1936 :lines 90 :wave 0 :status :done
   :callees ("N_BasisFindGlobalMax" "N_BasisIEval" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots"))
 (:c "N_BasisReparam" :lisp "basis-reparam" :file "NL_BasisBasic.cpp" :line 2202 :lines 18 :wave 0 :status :done
   :callees ("N_IntervalGetData" "N_KnotVectorGetKnots"))
 (:c "N_GetCompatibleKnotArray" :lisp "get-compatible-knot-array" :file "NL_BasisBasic.cpp" :line 2268 :lines 8 :wave 0 :status :done
   :callees ("N_GetCompatibleKnotVectorToTol"))
 (:c "N_GetCompatibleKnotVectorToTol" :lisp "get-compatible-knot-vector-to-tol" :file "NL_BasisBasic.cpp" :line 2483 :lines 116 :wave 0 :status :done
   :callees ("N_Alloc1dArrayKnotVectors" "N_AllocInt1dArray" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SetKnotIndex" "ST_notdone"))
 (:c "ST_notdone" :lisp "basis-basic--st-notdone" :file "NL_BasisBasic.cpp" :line 2695 :lines 19 :wave 0 :status :done :static t
   :callees ("N_KnotVectorGetKnots"))
 (:c "N_BasisFindIndexNodeArray" :lisp "basis-find-index-node-array" :file "NL_BasisBasic.cpp" :line 2757 :lines 20 :wave 0 :status :done
   :callees ("N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots"))
 (:c "N_BasisEval" :lisp "basis-eval" :file "NL_BasisBasic.cpp" :line 2838 :lines 49 :wave 0 :status :done
   :callees ("N_BasisFindSpan"))
 (:c "N_BasisDerivs" :lisp "basis-derivs" :file "NL_BasisBasic.cpp" :line 2987 :lines 105 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisFindSpan" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SwapIntegers"))
 (:c "N_BasisIEval" :lisp "basis-i-eval" :file "NL_BasisBasic.cpp" :line 3206 :lines 95 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_KnotVectorIsParamOutOfBounds"))
 (:c "N_BasisIDerivs" :lisp "basis-i-derivs" :file "NL_BasisBasic.cpp" :line 3394 :lines 158 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_AllocReal2dArray" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_KnotVectorIsParamOutOfBounds"))
 (:c "N_MakeKnotsCompatible" :lisp "make-knots-compatible" :file "NL_BasisBasic.cpp" :line 3665 :lines 57 :wave 0 :status :done
   :callees ("N_BasisReparam" "N_CreateInterval" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_BezGetDegreeElevationMatrix" :lisp "bez-get-degree-elevation-matrix" :file "NL_Bezier.cpp" :line 58 :lines 45 :wave 0 :status :done
   :callees ("N_CheckMemRealMatrix" "N_GetRealMatrixPtr"))
 (:c "N_BezFuncMultiplyBezCrvMatrix" :lisp "bez-func-multiply-bez-crv-matrix" :file "NL_Bezier.cpp" :line 184 :lines 36 :wave 0 :status :done
   :callees ("N_CheckMemRealMatrix" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs"))
 (:c "N_BezElevateDegree" :lisp "bez-elevate-degree" :file "NL_Bezier.cpp" :line 301 :lines 34 :wave 0 :status :done
   :callees ("N_BezGetDegreeElevationMatrix" "N_CopyCPt" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixIsNULL" "N_VectorBlendCPt"))
 (:c "N_BezFuncDegreeElevate" :lisp "bez-func-degree-elevate" :file "NL_Bezier.cpp" :line 402 :lines 33 :wave 0 :status :done
   :callees ("N_BezGetDegreeElevationMatrix" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixIsNULL"))
 (:c "N_BezFuncMultiplyBezFunc" :lisp "bez-func-multiply-bez-func" :file "NL_Bezier.cpp" :line 510 :lines 31 :wave 0 :status :done
   :callees ("N_BezFuncMultiplyBezCrvMatrix" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixIsNULL"))
 (:c "N_BezFuncMultiplyBezCrv" :lisp "bez-func-multiply-bez-crv" :file "NL_Bezier.cpp" :line 744 :lines 33 :wave 0 :status :done
   :callees ("N_BezFuncMultiplyBezCrvMatrix" "N_CopyCPt" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixIsNULL" "N_ScaleCPtXYZ" "N_VectorBlendCPt"))
 (:c "N_BezSrfElevateDegree" :lisp "bez-srf-elevate-degree" :file "NL_Bezier.cpp" :line 1509 :lines 60 :wave 0 :status :done
   :callees ("N_BezGetDegreeElevationMatrix" "N_CopyCPt" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixIsNULL" "N_VectorBlendCPt"))
 (:c "N_BezEvalOneBasis" :lisp "bez-eval-one-basis" :file "NL_Bezier.cpp" :line 1925 :lines 47 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_BezEvalBasis" :lisp "bez-eval-basis" :file "NL_Bezier.cpp" :line 2040 :lines 38 :wave 0 :status :port
   :callees ())
 (:c "N_BezDegreeReduceCoefs" :lisp "bez-degree-reduce-coefs" :file "NL_Bezier.cpp" :line 2476 :lines 15 :wave 0 :status :done
   :callees ())
 (:c "N_BezReduceDegree" :lisp "bez-reduce-degree" :file "NL_Bezier.cpp" :line 2540 :lines 58 :wave 0 :status :done
   :callees ("N_BezEvalOneBasis" "N_Combine2CPts" "N_CopyCPt" "N_DistCptCptHomo" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_BezConicArcFromPtsAndTangents" :lisp "bez-conic-arc-from-pts-and-tangents" :file "NL_Bezier.cpp" :line 2939 :lines 68 :wave 0 :status :done
   :callees ("N_CreateLinePtPt" "N_CreateLineStartDirVector" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_IsectLineLine" "N_VectorDiff" "N_VectorDot" "N_VectorScale"))
 (:c "N_BezSetConicWeight" :lisp "bez-set-conic-weight" :file "NL_Bezier.cpp" :line 3291 :lines 17 :wave 0 :status :done
   :callees ("N_DistPtPt"))
 (:c "N_BezConicSplit" :lisp "bez-conic-split" :file "NL_Bezier.cpp" :line 3350 :lines 29 :wave 0 :status :done
   :callees ("N_Combine2Pts" "N_CopyPt" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_BezReparam" :lisp "bez-reparam" :file "NL_Bezier.cpp" :line 4404 :lines 89 :wave 0 :status :done
   :callees ("N_AllocCPt2dArray" "N_AllocInt2dArray" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_InitNurbs" "N_PascalTriRow" "N_ScaleCPt" "N_VectorBlendCPt"))
 (:c "N_BezSplit" :lisp "bez-split" :file "NL_Bezier.cpp" :line 4568 :lines 50 :wave 0 :status :done
   :callees ("N_AllocCPt1dArray" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_BezCalcArcLength" :lisp "bez-calc-arc-length" :file "NL_Bezier.cpp" :line 4676 :lines 8 :wave 0 :status :done
   :callees ("ST_BezCalcArcLength"))
 (:c "ST_BezCalcArcLength" :lisp "bezier--st-bez-calc-arc-length" :file "NL_Bezier.cpp" :line 4713 :lines 45 :wave 0 :status :done :static t
   :callees ("N_AllocCPt1dArray" "N_BezSplit" "N_DistCPolygon" "N_DistCptCpt" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_BezCalcLength" :lisp "bez-calc-length" :file "NL_Bezier.cpp" :line 4832 :lines 23 :wave 0 :status :done
   :callees ("N_BezCalcArcLength" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_BezInterp4Pts" :lisp "bez-interp-4-pts" :file "NL_Bezier.cpp" :line 5173 :lines 53 :wave 0 :status :port
   :callees ("N_BezEvalBasis" "N_DistPtPt" "N_EndNurbs" "N_FloatOpIsBad" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_SetRealMatrix"))
 (:c "N_BezCenterRadiusFrom3CPts" :lisp "bez-center-radius-from-3-c-pts" :file "NL_Bezier.cpp" :line 5781 :lines 61 :wave 0 :status :done
   :callees ("N_CPtToPtEuclid" "N_Combine2Pts" "N_CopyPt" "N_DistPtLineSeg" "N_DistPtPt" "N_FloatOpIsBad" "N_PtFromXYZ" "N_VectorCreate" "N_VectorDiff" "N_VectorDirectedAngle" "N_VectorDot" "N_VectorNormalizeRef" "N_VectorPerpendicular" "N_VectorPtAlongVector" "N_VectorReverseInPlace"))
 (:c "N_CrvDerivsAtKnot" :lisp "crv-derivs-at-knot" :file "NL_CrvAdv.cpp" :line 181 :lines 34 :wave 3 :status :done
   :callees ("N_AllocReal2dArray" "N_BasisDerivs" "N_CopyCPt" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorIsParamOutOfBounds" "N_VectorBlendCPt"))
 (:c "N_CrvExtendByDist" :lisp "crv-extend-by-dist" :file "NL_CrvAdv.cpp" :line 630 :lines 4 :wave 3 :status :done :root t
   :callees ("ST_curxcdr"))
 (:c "ST_curxcdr" :lisp "crv-adv--st-curxcdr" :file "NL_CrvAdv.cpp" :line 640 :lines 580 :wave 3 :status :done :static t
   :callees ("N_AllocCPt1dArray" "N_AllocCrvArrays" "N_BasisFindSpanAndMult" "N_CPtFromWxWyWz" "N_CPtGetW" "N_CPtToPtEuclid" "N_CPtToWxWyWz" "N_Combine2CPts" "N_CopyCPt" "N_CptReflect" "N_CreatePlanePtNormal" "N_CrvArcLength" "N_CrvCopy" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvInitArrays" "N_CrvInsertKnot" "N_CrvRemoveKnot" "N_CrvSizeArrays" "N_EndNurbs" "N_FloatOpIsBad" "N_FreeCrv" "N_InitNurbs" "N_IsCrvRat" "N_PtToXYZ" "N_ScaleCPt" "N_VectorBlendCPt" "N_VectorBlendPt" "N_VectorDir" "N_VectorNormalize"))
 (:c "N_CrvExtendToPt" :lisp "crv-extend-to-pt" :file "NL_CrvAdv.cpp" :line 1656 :lines 544 :wave 3 :status :done :fenced t
   :callees ("N_AllocCrvArrays" "N_CPtGetW" "N_CPtSetW" "N_CPtSetZ" "N_CPtToPtEuclid" "N_CPtToWxWyWz" "N_Combine2CPts" "N_Combine2Pts" "N_CopyCPt" "N_CopyPt" "N_CreateLineStartDirVector" "N_CrvCopy" "N_CrvDerivs" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvInitArrays" "N_CrvInsertKnot" "N_CrvIs3d" "N_CrvRemoveKnot" "N_CrvSizeArrays" "N_DistCptCpt" "N_DistPtPt" "N_EndNurbs" "N_FloatOpIsBad" "N_FreeCrv" "N_InitNurbs" "N_IsCrvRat" "N_IsectLineLine" "N_ProjectPtLine" "N_PtToCPt" "N_ScaleCPt" "N_ScalePt" "N_VectorBlendCPt" "N_VectorCross" "N_VectorDir" "N_VectorDot" "N_VectorMagnitude" "N_VectorNormalize"))
 (:c "N_CrvGetMax2ndDeriv" :lisp "crv-get-max-2nd-deriv" :file "NL_CrvAdv.cpp" :line 3059 :lines 48 :wave 3 :status :done :root t
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_CrvDecomposeBez" "N_CrvDerivs" "N_CrvGetDegree" "N_CrvGetKnots" "N_EndNurbs" "N_InitNurbs" "N_PtMagnitude"))
 (:c "N_CrvPlanarOffsetGetDeriv" :lisp "crv-planar-offset-get-deriv" :file "NL_CrvAdv.cpp" :line 3184 :lines 35 :wave 3 :status :done
   :callees ("N_CrvDerivs" "N_CrvGetCurvatureDeriv" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_VectorCombine" "N_VectorCross" "N_VectorDot" "N_VectorMagnitude" "N_VectorScale"))
 (:c "N_CrvOffsetGetMax2ndDeriv" :lisp "crv-offset-get-max-2nd-deriv" :file "NL_CrvAdv.cpp" :line 3279 :lines 46 :wave 3 :status :done
   :callees ("N_AllocReal1dArray" "N_CrvDecomposeBez" "N_CrvGetDegree" "N_CrvGetKnots" "N_CrvPlanarOffsetGetDeriv" "N_EndNurbs" "N_InitNurbs" "N_PtMagnitude"))
 (:c "N_CrvRemoveAllKnotsArcLen" :lisp "crv-remove-all-knots-arc-len" :file "NL_CrvAdv.cpp" :line 4807 :lines 27 :wave 3 :status :done :root t :fenced t
   :callees ("N_CrvRemoveKnots" "N_EndNurbs" "N_InitNurbs" "N_SrfReparamMultKnots"))
 (:c "N_CrvRemoveAllKnots" :lisp "crv-remove-all-knots" :file "NL_CrvAdv.cpp" :line 4901 :lines 228 :wave 3 :status :done :root t
   :callees ("N_AllocCPt1dArray" "N_AllocInt1dArray" "N_AllocReal1dArray" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_CrvCompress" "N_CrvCopy" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetCPtsKnotVectorAndKnots" "N_CrvGetMinMaxWeightsAndPts" "N_CrvRemoveKnotMaxErr" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat"))
 (:c "N_CrvsMakeCompatibleApprox" :lisp "crvs-make-compatible-approx" :file "NL_CrvAdv.cpp" :line 5733 :lines 357 :wave 3 :status :done :root t
   :callees ("N_AllocCPt2dArray" "N_AllocInt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisFindAllSpanMaxima" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_CrvCompress" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvGetMinMaxWeightsAndPts" "N_CrvRemoveKnotMaxErr" "N_CrvSetSizeIndices" "N_CrvsMakeCompatibleKnotTol" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat"))
 (:c "N_CrvRemoveKnotsDerivConstraints" :lisp "crv-remove-knots-deriv-constraints" :file "NL_CrvAdv.cpp" :line 8121 :lines 295 :wave 3 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocInt1dArray" "N_AllocReal1dArray" "N_BasisFindAllSpanMaxima" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_CrvCompress" "N_CrvCopy" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetCPtsKnotVectorAndKnots" "N_CrvGetMinMaxWeightsAndPts" "N_CrvRemoveKnotMaxErr" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat"))
 (:c "N_CrvMakeCompatibleWithSrf" :lisp "crv-make-compatible-with-srf" :file "NL_CrvAdv.cpp" :line 8603 :lines 72 :wave 3 :status :done
   :callees ("N_BasisReparam" "N_CreateInterval" "N_CrvElevateDegree" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvNonRatToRat" "N_CrvRefine" "N_EndNurbs" "N_GetCompatibleKnotArray" "N_InitNurbs" "N_IsCrvRat" "N_IsSrfRat" "N_KnotVectorGetKnots" "N_SrfElevateDegree" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_SrfInsertKnots" "N_SrfNonRatToRat"))
 (:c "N_CrvGetRatSecondDeriv" :lisp "crv-get-rat-second-deriv" :file "NL_CrvAdv.cpp" :line 9134 :lines 143 :wave 3 :status :done :root t :fenced t
   :callees ("N_AllocKnotVectorAndArray" "N_BasisGetSpanCount" "N_BasisIncreaseKnotMult" "N_CFuncGetDegree" "N_CFuncGetKnotVector" "N_CFuncInitArrays" "N_CreateCrvFromNumAndDenom" "N_CrvElevateDegree" "N_CrvFuncDeriv" "N_CrvFuncMultiplyConstant" "N_CrvFuncMultiplyCrv" "N_CrvFuncMultiplyCrvFunc" "N_CrvFuncRefine" "N_CrvFuncSumDiffCrvFunc" "N_CrvGetDegree" "N_CrvGetKnotVector" "N_CrvGetNumAndDenom" "N_CrvInitArrays" "N_CrvNonRatGetFirstDeriv" "N_CrvRefine" "N_CrvSumDiffCrv" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat" "N_KnotVectorGetKnots"))
 (:c "N_CrvSumDiffCrv" :lisp "crv-sum-diff-crv" :file "NL_CrvAdv.cpp" :line 10058 :lines 155 :wave 3 :status :done
   :callees ("N_AllocKnotVectorAndArray" "N_BasisGetSpanCount" "N_CFuncInitArrays" "N_CreateCrvFromNumAndDenom" "N_CrvCopy" "N_CrvElevateDegree" "N_CrvFuncMultiplyCrv" "N_CrvFuncMultiplyCrvFunc" "N_CrvGetCPts" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvGetNumAndDenom" "N_CrvInitArrays" "N_CrvNonRatToRat" "N_CrvRatToNonRat" "N_CrvRefine" "N_CrvSizeArrays" "N_Diff2CPts" "N_EndNurbs" "N_GetCompatibleKnotArrayMult" "N_InitNurbs" "N_IsCrvRat" "N_KnotVectorGetKnots" "N_Sum2CPts"))
 (:c "N_CrvNonRatGetFirstDeriv" :lisp "crv-non-rat-get-first-deriv" :file "NL_CrvAdv.cpp" :line 10995 :lines 52 :wave 3 :status :done :root t
   :callees ("N_AllocCPt1dArray" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvSizeArrays" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_IsCrvRat"))
 (:c "N_CrvRatGetFirstDeriv" :lisp "crv-rat-get-first-deriv" :file "NL_CrvAdv.cpp" :line 11128 :lines 73 :wave 3 :status :done :root t
   :callees ("N_AllocKnotVectorAndArray" "N_BasisGetSpanCount" "N_BasisIncreaseKnotMult" "N_CFuncGetDegree" "N_CFuncGetKnotVector" "N_CFuncInitArrays" "N_CreateCrvFromNumAndDenom" "N_CrvElevateDegree" "N_CrvFuncDeriv" "N_CrvFuncMultiplyCrv" "N_CrvFuncMultiplyCrvFunc" "N_CrvFuncRefine" "N_CrvGetNumAndDenom" "N_CrvInitArrays" "N_CrvNonRatGetFirstDeriv" "N_CrvSumDiffCrv" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat" "N_KnotVectorGetKnots"))
 (:c "N_SrfReparamFunc" :lisp "srf-reparam-func" :file "NL_CrvAdv.cpp" :line 13176 :lines 321 :wave 3 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocCrvArrays" "N_AllocKnotVectorAndArray" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_BezReparam" "N_CFuncGetData" "N_CFuncGetDegree" "N_CFuncGetKnotVector" "N_CFuncInitArrays" "N_Combine2CPts" "N_CopyCPt" "N_CrvCompress" "N_CrvFuncRefine" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetDegree" "N_CrvGetKnotVector" "N_CrvInitArrays" "N_CrvRefine" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_EndNurbs" "N_FreeCrv" "N_InitNurbs" "N_KnotVectorGetKnots" "N_MapKnotsBetweenCrvFuncAndKnotVector"))
 (:c "N_SrfReparamMultKnots" :lisp "srf-reparam-mult-knots" :file "NL_CrvAdv.cpp" :line 13712 :lines 79 :wave 3 :status :done :root t
   :callees ("N_CopyCPt" "N_CrvArcLength" "N_CrvCopy" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvInitArrays" "N_CrvSizeArrays" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_ApproxContinuousCrvWithArcs" :lisp "approx-continuous-crv-with-arcs" :file "NL_CrvApprox.cpp" :line 239 :lines 252 :wave 5 :status :done :root t
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_ApproxCrvWithPolyline" "N_CPtToPtEuclid" "N_Combine2Pts" "N_CopyCPt" "N_CrvAreArraysNULL" "N_CrvCompress" "N_CrvDerivs" "N_CrvExpand" "N_CrvGetArraySizes" "N_CrvGetCPtsAndKnots" "N_CrvGetDegree" "N_CrvIsPolygonWithinTol" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_DistPtPt" "N_EndNurbs" "N_FitArcToEndPtsAndTangents" "N_FitCalcCrvParamValues" "N_InitNurbs" "N_Realloc1dRealArray" "N_VectorCopy" "N_Weight"))
 (:c "N_CrvIsPolygonWithinTol" :lisp "crv-is-polygon-within-tol" :file "NL_CrvApprox.cpp" :line 870 :lines 86 :wave 5 :status :done :root t
   :callees ("N_AllocInt1dArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_BezCenterRadiusFrom3CPts" "N_DistMaxLineArc" "N_DistPtArc" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_ApproxNurbsWithNonRatCrv" :lisp "approx-nurbs-with-non-rat-crv" :file "NL_CrvApprox.cpp" :line 2698 :lines 161 :wave 5 :status :done :root t
   :callees ("N_AllocInt1dArray" "N_AllocKnotVectorAndArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_CreateInterval" "N_Crv3dTo2d" "N_CrvArcLength" "N_CrvCompress" "N_CrvDecomposeBez" "N_CrvDerivs" "N_CrvEval" "N_CrvGetArraySizes" "N_CrvGetDegree" "N_CrvGetKnots" "N_CrvGetMax2ndDeriv" "N_CrvInitArrays" "N_CrvIs3d" "N_CrvRemoveKnotsDerivConstraints" "N_CrvReparamToInterval" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCalcKnotVectorEndDerivs" "N_FitCrvCalcKnotVector" "N_FitCrvInterpGivenParams" "N_FitCrvKnotsAndDerivs" "N_FloatOpIsBad" "N_InitNurbs" "N_VectorCopy" "N_VectorNormalizeRef" "N_VectorScale"))
 (:c "N_CrvIs3d" :lisp "crv-is-3d" :file "NL_CrvBasic.cpp" :line 38 :lines 15 :wave 0 :status :done
   :callees ("N_CPtToWxWyWz"))
 (:c "N_CrvIsSized" :lisp "crv-is-sized" :file "NL_CrvBasic.cpp" :line 432 :lines 12 :wave 0 :status :vanish
   :callees ("N_CrvGetArraySizes" "N_ErrSet"))
 (:c "N_AllocCrv" :lisp "alloc-crv" :file "NL_CrvBasic.cpp" :line 686 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocCrvAndArrays" :lisp "alloc-crv-and-arrays" :file "NL_CrvBasic.cpp" :line 761 :lines 17 :wave 0 :status :vanish
   :callees ("N_AllocCPolygonAndArray" "N_AllocCrv" "N_AllocKnotVectorAndArray" "N_CrvFromCPolygonAndKnotVector"))
 (:c "N_Alloc1dArrayCrvs" :lisp "alloc-1d-array-crvs" :file "NL_CrvBasic.cpp" :line 832 :lines 35 :wave 0 :status :vanish
   :callees ("N_AllocCrvAndArrays" "N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocArrayCrvPtrs" :lisp "alloc-array-crv-ptrs" :file "NL_CrvBasic.cpp" :line 918 :lines 24 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocArrayCrvPtrsAndData" :lisp "alloc-array-crv-ptrs-and-data" :file "NL_CrvBasic.cpp" :line 992 :lines 33 :wave 0 :status :vanish :root t
   :callees ("N_AllocCrv" "N_CrvInitArrays" "N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_CrvFromCPolygonAndKnotVector" :lisp "crv-from-c-polygon-and-knot-vector" :file "NL_CrvBasic.cpp" :line 1335 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_CrvFromCPtsAndKnots" :lisp "crv-from-c-pts-and-knots" :file "NL_CrvBasic.cpp" :line 1384 :lines 15 :wave 0 :status :done :root t
   :callees ("N_AllocCPolygon" "N_AllocKnotVector" "N_CPolygonFromCPts" "N_CrvFromCPolygonAndKnotVector" "N_KnotVectorFromRealArray"))
 (:c "N_AllocCrvArrays" :lisp "alloc-crv-arrays" :file "NL_CrvBasic.cpp" :line 1523 :lines 18 :wave 0 :status :vanish
   :callees ("N_AllocCPolygonAndArray" "N_AllocKnotVectorAndArray" "N_CrvFromCPolygonAndKnotVector"))
 (:c "N_CrvSizeArrays" :lisp "crv-size-arrays" :file "NL_CrvBasic.cpp" :line 1593 :lines 24 :wave 0 :status :vanish
   :callees ("N_AllocCrvArrays" "N_CrvAreArraysNULL" "N_CrvIsSized" "N_CrvSetSizeIndices"))
 (:c "N_CrvDetachPolygonKnot" :lisp "crv-detach-polygon-knot" :file "NL_CrvBasic.cpp" :line 1741 :lines 6 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvGetCPtsDegreeAndKnots" :lisp "crv-get-c-pts-degree-and-knots" :file "NL_CrvBasic.cpp" :line 1783 :lines 25 :wave 0 :status :vanish
   :callees ("N_CrvAreArraysNULL"))
 (:c "N_CrvGetArraySizes" :lisp "crv-get-array-sizes" :file "NL_CrvBasic.cpp" :line 1836 :lines 5 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvGetConicData" :lisp "crv-get-conic-data" :file "NL_CrvBasic.cpp" :line 1863 :lines 41 :wave 0 :status :done
   :callees ())
 (:c "N_CrvGetCPts" :lisp "crv-get-c-pts" :file "NL_CrvBasic.cpp" :line 1950 :lines 8 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvGetCPtsAndKnots" :lisp "crv-get-c-pts-and-knots" :file "NL_CrvBasic.cpp" :line 1988 :lines 8 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvGetCPtsKnotVectorAndKnots" :lisp "crv-get-c-pts-knot-vector-and-knots" :file "NL_CrvBasic.cpp" :line 2028 :lines 6 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvSetCPtsAndKnots" :lisp "crv-set-c-pts-and-knots" :file "NL_CrvBasic.cpp" :line 2063 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_CrvGetKnots" :lisp "crv-get-knots" :file "NL_CrvBasic.cpp" :line 2097 :lines 5 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_CrvGetParamBounds" :lisp "crv-get-param-bounds" :file "NL_CrvBasic.cpp" :line 2131 :lines 8 :wave 0 :status :done :root t
   :callees ())
 (:c "N_CrvGetKnotVector" :lisp "crv-get-knot-vector" :file "NL_CrvBasic.cpp" :line 2166 :lines 4 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvGetDegree" :lisp "crv-get-degree" :file "NL_CrvBasic.cpp" :line 2195 :lines 4 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_CrvCompress" :lisp "crv-compress" :file "NL_CrvBasic.cpp" :line 2235 :lines 26 :wave 0 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_CopyCPt" "N_CrvGetArraySizes" "N_CrvGetCPtsAndKnots" "N_CrvSetCPtsAndKnots" "N_FreeCPt1dArray" "N_FreeReal1dArray"))
 (:c "N_CrvExpand" :lisp "crv-expand" :file "NL_CrvBasic.cpp" :line 2316 :lines 31 :wave 0 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_CopyCPt" "N_CrvGetArraySizes" "N_CrvGetCPtsAndKnots" "N_CrvGetDegree" "N_CrvSetCPtsAndKnots" "N_CrvSetSizeIndices" "N_FreeCPt1dArray" "N_FreeReal1dArray"))
 (:c "N_CrvSetSizeIndices" :lisp "crv-set-size-indices" :file "NL_CrvBasic.cpp" :line 2437 :lines 6 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvInitArrays" :lisp "crv-init-arrays" :file "NL_CrvBasic.cpp" :line 2472 :lines 6 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_CrvAreArraysNULL" :lisp "crv-are-arrays-null" :file "NL_CrvBasic.cpp" :line 2506 :lines 15 :wave 0 :status :vanish
   :callees ("N_CrvDetachPolygonKnot"))
 (:c "N_IsCrvRat" :lisp "is-crv-rat" :file "NL_CrvBasic.cpp" :line 2554 :lines 16 :wave 0 :status :done :root t
   :callees ("N_CPtGetW" "N_CrvGetCPts"))
 (:c "N_CrvIsClosed" :lisp "crv-is-closed" :file "NL_CrvBasic.cpp" :line 2605 :lines 19 :wave 0 :status :done
   :callees ("N_CPtToPtEuclid" "N_CrvGetCPts" "N_DistPtPt"))
 (:c "N_CrvIsDegen" :lisp "crv-is-degen" :file "NL_CrvBasic.cpp" :line 2660 :lines 31 :wave 0 :status :done
   :callees ("N_CPtToPtEuclid" "N_CopyPt" "N_CrvGetCPts" "N_DistPtPt" "N_ScalePt" "N_Sum2Pts"))
 (:c "N_CrvIsLine" :lisp "crv-is-line" :file "NL_CrvBasic.cpp" :line 2736 :lines 50 :wave 0 :status :done
   :callees ("N_CPtToPtEuclid" "N_CrvGetCPts" "N_DistPerpPtLineSeg" "N_DistPtPt" "N_VectorDir" "N_VectorMagnitude" "N_VectorNormalizeRef" "N_VectorPtAlongVector" "N_VectorReverseInPlace"))
 (:c "N_CrvsAreCoincident" :lisp "crvs-are-coincident" :file "NL_CrvBasic.cpp" :line 2843 :lines 33 :wave 0 :status :done
   :callees ("N_AllocArrayCrvPtrs" "N_CrvGetCPts" "N_CrvsMakeCompatible" "N_DistCptCptHomo"))
 (:c "N_CrvsAreCombatible" :lisp "crvs-are-combatible" :file "NL_CrvBasic.cpp" :line 3029 :lines 78 :wave 0 :status :done
   :callees ("N_CrvGetArraySizes" "N_CrvGetKnots" "N_CrvIs3d" "N_IsCrvRat"))
 (:c "N_CrvCopy" :lisp "crv-copy" :file "NL_CrvBasic.cpp" :line 3178 :lines 19 :wave 0 :status :done
   :callees ("N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvSizeArrays"))
 (:c "N_CrvGetNumAndDenom" :lisp "crv-get-num-and-denom" :file "NL_CrvBasic.cpp" :line 3256 :lines 37 :wave 0 :status :done
   :callees ("N_CFuncSizeArrays" "N_CPtFromWxWyWz" "N_CPtToWxWyWz" "N_CrvFuncCntrlValKnots" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvSizeArrays" "N_IsCrvRat"))
 (:c "N_CreateCrvFromNumAndDenom" :lisp "create-crv-from-num-and-denom" :file "NL_CrvBasic.cpp" :line 3493 :lines 33 :wave 0 :status :done
   :callees ("N_CFuncGetData" "N_CPtFromWxWyWz" "N_CPtToWxWyWz" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvSizeArrays"))
 (:c "N_CrvGetEPolygon" :lisp "crv-get-e-polygon" :file "NL_CrvBasic.cpp" :line 4514 :lines 22 :wave 0 :status :done
   :callees ("N_AllocPt1dArray" "N_CPtToPtEuclid" "N_CrvGetCPts" "N_EPolygonFromPts" "N_ErrSet"))
 (:c "N_FreeCrv" :lisp "free-crv" :file "NL_CrvBasic.cpp" :line 4586 :lines 18 :wave 0 :status :vanish
   :callees ("N_CrvDetachPolygonKnot" "N_CrvGetCPtsAndKnots" "N_FreeCPolygon" "N_FreeCPt1dArray" "N_FreeKnotVector" "N_FreeReal1dArray"))
 (:c "N_Crv3dTo2d" :lisp "crv-3d-to-2d" :file "NL_CrvBasic.cpp" :line 4717 :lines 10 :wave 0 :status :done
   :callees ("N_CPtSetZ" "N_CrvGetCPts"))
 (:c "N_Crv2dTo3d" :lisp "crv-2d-to-3d" :file "NL_CrvBasic.cpp" :line 4765 :lines 10 :wave 0 :status :done
   :callees ("N_CPtSetZ" "N_CrvGetCPts"))
 (:c "N_CrvRatToNonRat" :lisp "crv-rat-to-non-rat" :file "NL_CrvBasic.cpp" :line 4813 :lines 10 :wave 0 :status :done
   :callees ("N_CPtSetW" "N_CrvGetCPts"))
 (:c "N_CrvNonRatToRat" :lisp "crv-non-rat-to-rat" :file "NL_CrvBasic.cpp" :line 4860 :lines 10 :wave 0 :status :done
   :callees ("N_CPtSetW" "N_CrvGetCPts"))
 (:c "N_CrvsMakeRatCompatible" :lisp "crvs-make-rat-compatible" :file "NL_CrvBasic.cpp" :line 4910 :lines 45 :wave 0 :status :done
   :callees ("N_Crv2dTo3d" "N_CrvIs3d" "N_CrvNonRatToRat" "N_IsCrvRat"))
 (:c "N_CrvReparamToInterval" :lisp "crv-reparam-to-interval" :file "NL_CrvBasic.cpp" :line 5005 :lines 20 :wave 0 :status :done
   :callees ("N_CrvGetDegree" "N_CrvGetKnots" "N_IntervalGetData"))
 (:c "N_CrvIsPlanar" :lisp "crv-is-planar" :file "NL_CrvBasic.cpp" :line 5364 :lines 94 :wave 0 :status :done
   :callees ("N_AllocPt1dArray" "N_CPtToPtEuclid" "N_CopyPt" "N_CrvGetCPts" "N_CrvIsDegen" "N_CrvIsLine" "N_DistSignedPtPlane" "N_FreePt1dArray" "N_PlanePtNormalToImplicit" "N_ScalePt" "N_Sum2Pts" "N_VectorCross" "N_VectorDir" "N_VectorDot" "N_VectorNormalize" "N_VectorScale" "N_VectorSum" "N_VectorsAngle"))
 (:c "N_CreateQuadraticArcTol" :lisp "create-quadratic-arc-tol" :file "NL_CrvCommon.cpp" :line 75 :lines 104 :wave 3 :status :done
   :callees ("N_Combine2Pts" "N_CopyCPt" "N_CopyPt" "N_CreateLineStartDirVector" "N_CrvGetCPtsAndKnots" "N_CrvGetConicData" "N_CrvSizeArrays" "N_DistPtPt" "N_IsectLineLine" "N_TranslateSum2Pts" "N_VectorNormalizeRef" "N_Weight"))
 (:c "N_CreateQuarticArcTol" :lisp "create-quartic-arc-tol" :file "NL_CrvCommon.cpp" :line 311 :lines 126 :wave 3 :status :done
   :callees ("N_Combine2CPts" "N_Combine2Pts" "N_Combine4CPts" "N_CopyCPt" "N_CopyPt" "N_CreateLineStartDirVector" "N_CrvGetCPtsAndKnots" "N_CrvGetConicData" "N_CrvSizeArrays" "N_IsectLineLine" "N_TranslateSum2Pts" "N_VectorNormalizeRef" "N_Weight"))
 (:c "N_CreateQuinticArcTol" :lisp "create-quintic-arc-tol" :file "NL_CrvCommon.cpp" :line 581 :lines 87 :wave 3 :status :done
   :callees ("N_Combine2Pts" "N_CopyCPt" "N_CrvClosestPt" "N_CrvExtractCrvSeg" "N_CrvGetCPtsAndKnots" "N_CrvGetConicData" "N_CrvSizeArrays" "N_TranslateSum2Pts" "N_VectorNormalizeRef" "N_Weight"))
 (:c "N_CreateCircArc" :lisp "create-circ-arc" :file "NL_CrvCommon.cpp" :line 745 :lines 16 :wave 3 :status :done :root t
   :callees ("N_CreateCircArcTol"))
 (:c "N_CreateCircArcTol" :lisp "create-circ-arc-tol" :file "NL_CrvCommon.cpp" :line 823 :lines 37 :wave 3 :status :done
   :callees ("N_CreateQuadraticArcTol" "N_CreateQuarticArcTol" "N_CreateQuinticArcTol"))
 (:c "N_CalcCircWeights" :lisp "calc-circ-weights" :file "NL_CrvCommon.cpp" :line 919 :lines 22 :wave 3 :status :done
   :callees ("N_CPtGetW" "N_CreateCircArc" "N_CrvGetCPts" "N_CrvInitArrays" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_CreateEllipticalArc" :lisp "create-elliptical-arc" :file "NL_CrvCommon.cpp" :line 971 :lines 8 :wave 3 :status :done :root t
   :callees ("N_CreateEllipticalArcAngTypeTol"))
 (:c "N_CreateEllipticalArcAngTypeTol" :lisp "create-elliptical-arc-ang-type-tol" :file "NL_CrvCommon.cpp" :line 1090 :lines 77 :wave 3 :status :done
   :callees ("N_CPtGetW" "N_CPtToPtEuclid" "N_Combine2Pts" "N_CreateCircArcTol" "N_CreateLineStartDirVector" "N_CrvGetCPts" "N_ProjectPtLine" "N_VectorNormalizeRef" "N_Weight"))
 (:c "N_CreateConicArc" :lisp "create-conic-arc" :file "NL_CrvCommon.cpp" :line 1269 :lines 96 :wave 3 :status :done :root t
   :callees ("N_BezConicArcFromPtsAndTangents" "N_BezConicSplit" "N_CopyPt" "N_CrvGetCPtsAndKnots" "N_CrvSizeArrays" "N_VectorDot" "N_VectorMagnitude" "N_Weight"))
 (:c "N_FitArcToPts" :lisp "fit-arc-to-pts" :file "NL_CrvFit.cpp" :line 89 :lines 140 :wave 4 :status :done :root t :fenced t
   :callees ("N_AllocReal1dArray" "N_CreateCircArc" "N_CreateCircFrom2dPts" "N_CrvLineFromPtAndVector" "N_EndNurbs" "N_InitNurbs" "N_LineFitPts" "N_ProjectPtLineParam" "N_ShellSortReal" "N_VectorCreate" "N_VectorDiff" "N_VectorDirectedAngle" "N_VectorSum"))
 (:c "N_FitArcToEndPtsAndTangents" :lisp "fit-arc-to-end-pts-and-tangents" :file "NL_CrvFit.cpp" :line 353 :lines 200 :wave 4 :status :done
   :callees ("N_Combine2Pts" "N_DistPtLineSeg" "N_DistPtPt" "N_PtToXYZ" "N_VectorCross" "N_VectorDiff" "N_VectorDot" "N_VectorNormalize" "N_VectorNormalizeRef" "N_VectorPerpendicular" "N_VectorPtAlongVector" "N_VectorReverse" "N_VectorSum" "N_Weight" "ST_CalcCPtsCircle"))
 (:c "N_FitCalcCrvParamValues" :lisp "fit-calc-crv-param-values" :file "NL_CrvFit.cpp" :line 1361 :lines 81 :wave 4 :status :done
   :callees ("N_AllocReal1dArray" "N_DistCptCptHomo" "N_DistPtPt" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_FitFuncInterpGivenParams" :lisp "fit-func-interp-given-params" :file "NL_CrvFit.cpp" :line 1845 :lines 58 :wave 4 :status :done
   :callees ("N_AllocReal1dArray" "N_BasisEval" "N_CFuncSizeArrays" "N_CrvFuncCntrlValKnotVector" "N_EndNurbs" "N_FitCrvCalcKnotVector" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixLuDecompose" "N_RealMatrixRightForBack" "N_SetRealMatrix"))
 (:c "N_FitCrvCalcKnotVector" :lisp "fit-crv-calc-knot-vector" :file "NL_CrvFit.cpp" :line 2000 :lines 26 :wave 4 :status :done
   :callees ("N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_FitCrvInterp" :lisp "fit-crv-interp" :file "NL_CrvFit.cpp" :line 2092 :lines 66 :wave 4 :status :done :root t
   :callees ("N_AllocReal1dArray" "N_BasisEval" "N_CrvGetCPtsKnotVectorAndKnots" "N_CrvSizeArrays" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCrvCalcKnotVector" "N_GetRealMatrixPtr" "N_InitNurbs" "N_PtToCPt" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_SetRealMatrix"))
 (:c "N_FitCrvInterpGivenParams" :lisp "fit-crv-interp-given-params" :file "NL_CrvFit.cpp" :line 2273 :lines 86 :wave 4 :status :done
   :callees ("N_AllocReal1dArray" "N_BasisEval" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvSizeArrays" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_KnotVectorGetKnots" "N_PtToCPt" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_SetRealMatrix"))
 (:c "N_FitCalcKnotVectorEndDerivs" :lisp "fit-calc-knot-vector-end-derivs" :file "NL_CrvFit.cpp" :line 2458 :lines 41 :wave 4 :status :done
   :callees ("N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_FitCalcKnotVectorDeriv" :lisp "fit-calc-knot-vector-deriv" :file "NL_CrvFit.cpp" :line 2562 :lines 70 :wave 4 :status :done
   :callees ("N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_FitCrvDerivs" :lisp "fit-crv-derivs" :file "NL_CrvFit.cpp" :line 2716 :lines 110 :wave 4 :status :done :root t
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_BasisEval" "N_CopyPt" "N_CrvDerivs" "N_CrvGetCPtsKnotVectorAndKnots" "N_CrvInitArrays" "N_CrvSizeArrays" "N_DistPolygon" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCalcKnotVectorEndDerivs" "N_FitCrvInterp" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_ScalePt" "N_SetRealMatrix" "N_VectorCopy" "N_VectorNormalizeRef" "N_VectorScale"))
 (:c "N_FitCrvKnotsAndDerivs" :lisp "fit-crv-knots-and-derivs" :file "NL_CrvFit.cpp" :line 3063 :lines 116 :wave 4 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_BasisEval" "N_CopyCPt" "N_CopyPt" "N_CrvGetCPtsAndKnots" "N_CrvSizeArrays" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_KnotVectorGetKnots" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_ScaleCPt" "N_ScalePt" "N_SetRealMatrix"))
 (:c "N_FitCrvKnotsAndDeriv" :lisp "fit-crv-knots-and-deriv" :file "NL_CrvFit.cpp" :line 3313 :lines 165 :wave 4 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_BasisEval" "N_CopyCPt" "N_CopyPt" "N_CrvGetCPtsAndKnots" "N_CrvSizeArrays" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_KnotVectorGetKnots" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_ScaleCPt" "N_ScalePt" "N_SetRealMatrix"))
 (:c "N_FitCubicSplineInterp" :lisp "fit-cubic-spline-interp" :file "NL_CrvFit.cpp" :line 3615 :lines 127 :wave 4 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_BasisEval" "N_Combine2CPts" "N_Combine2Pts" "N_CrvDerivs" "N_CrvGetCPtsKnotVectorAndKnots" "N_CrvInitArrays" "N_CrvSizeArrays" "N_DistPolygon" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCrvInterp" "N_InitNurbs" "N_PtToCPt" "N_VectorBlendCPt" "N_VectorCopy" "N_VectorNormalizeRef" "N_VectorScale"))
 (:c "N_FitCrvCPtsFromSrfData" :lisp "fit-crv-c-pts-from-srf-data" :file "NL_CrvFit.cpp" :line 3877 :lines 65 :wave 4 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocPt1dArray" "N_CopyCPt" "N_CopyPt" "N_EndNurbs" "N_InitNurbs" "N_RealMatrixForBack"))
 (:c "N_FitCalcKnotVectorCrvApprox" :lisp "fit-calc-knot-vector-crv-approx" :file "NL_CrvFit.cpp" :line 4023 :lines 57 :wave 4 :status :done
   :callees ("N_AllocReal1dArray" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_FitCrvApprox" :lisp "fit-crv-approx" :file "NL_CrvFit.cpp" :line 4482 :lines 206 :wave 4 :status :done :root t
   :callees ("N_1dPtSetIsClosed" "N_AllocCrvArrays" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_CopyPt" "N_CrvAreArraysNULL" "N_CrvClosestPt" "N_CrvCompress" "N_CrvDerivs" "N_CrvGetArraySizes" "N_CrvGetKnotVector" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_DistPtPt" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCalcKnotVectorEndDerivs" "N_FitCrvApproxKnots" "N_FitCrvCalcKnotVector" "N_FitCrvInterpGivenParams" "N_FitCrvKnotsAndDerivs" "N_FitCrvWeightedLstSqKnots" "N_FitRemoveKnots" "N_InitNurbs" "N_KnotVectorGetKnots" "N_Pts1dGetMaxExtent" "N_VectorCopy"))
 (:c "N_FitCrvApproxKnots" :lisp "fit-crv-approx-knots" :file "NL_CrvFit.cpp" :line 5461 :lines 203 :wave 4 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocInt1dArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisEval" "N_BasisIEval" "N_CopyCPt" "N_CopyPt" "N_CrvGetCPtsAndKnots" "N_CrvSizeArrays" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_KnotVectorGetKnots" "N_PtToCPt" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_SetRealMatrix" "N_TranslateSum2CPts" "N_TranslateSum2Pts" "N_VectorBlendCPt" "N_VectorBlendPt"))
 (:c "N_FitCrvFirstDeriv" :lisp "fit-crv-first-deriv" :file "NL_CrvFit.cpp" :line 7689 :lines 101 :wave 4 :status :done :root t
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisDerivs" "N_CopyPt" "N_CrvGetCPtsAndKnots" "N_CrvGetKnotVector" "N_CrvSizeArrays" "N_DistPolygon" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCalcKnotsDerivs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_ScalePt" "N_SetRealMatrix" "N_VectorNormalizeRef" "N_VectorScale"))
 (:c "N_FitCrvFirstDerivAndKnots" :lisp "fit-crv-first-deriv-and-knots" :file "NL_CrvFit.cpp" :line 7923 :lines 121 :wave 4 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocPt1dArray" "N_AllocReal2dArray" "N_BasisDerivs" "N_CopyCPt" "N_CopyPt" "N_CrvGetCPtsAndKnots" "N_CrvSizeArrays" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_KnotVectorGetKnots" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_ScaleCPt" "N_ScalePt" "N_SetRealMatrix"))
 (:c "N_FitCrvKnotsAndTangents" :lisp "fit-crv-knots-and-tangents" :file "NL_CrvFit.cpp" :line 8207 :lines 251 :wave 4 :status :done
   :callees ("N_AllocKnotVectorAndArray" "N_AllocReal1dArray" "N_DistPolygon" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCalcKnotVectorDeriv" "N_FitCalcKnotVectorEndDerivs" "N_FitCrvCalcKnotVector" "N_FitCrvInterpGivenParams" "N_FitCrvKnotsAndDeriv" "N_FitCrvKnotsAndDerivs" "N_FreeKnotVector" "N_InitNurbs" "N_KnotVectorGetKnots" "N_KnotVectorIsValid" "N_VectorNormalizeRef" "N_VectorScale"))
 (:c "N_FitCrvConics" :lisp "fit-crv-conics" :file "NL_CrvFit.cpp" :line 10293 :lines 309 :wave 4 :status :done :root t :fenced t
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_BezSetConicWeight" "N_Combine2Pts" "N_CopyPt" "N_CreateLineStartDirVector" "N_CrvAreArraysNULL" "N_CrvCompress" "N_CrvGetCPtsAndKnots" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_DistPtPt" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FloatOpIsBad" "N_InitNurbs" "N_IsectLineLine" "N_VectorCombine" "N_VectorCopy" "N_VectorCross" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude" "N_VectorNormalize" "N_VectorNormalizeRef" "N_VectorPtAlongVector" "N_Weight"))
 (:c "N_FitRemoveKnots" :lisp "fit-remove-knots" :file "NL_CrvFit.cpp" :line 11049 :lines 215 :wave 4 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocInt1dArray" "N_AllocReal1dArray" "N_BasisIEval" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvRemoveKnotMaxErr" "N_CrvSetSizeIndices" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_FitCrvWeightedLstSqKnots" :lisp "fit-crv-weighted-lst-sq-knots" :file "NL_CrvFit.cpp" :line 13320 :lines 445 :wave 4 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocInt1dArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisDerivs" "N_CopyCPt" "N_CopyPt" "N_CrvGetCPtsAndKnots" "N_CrvSizeArrays" "N_Diff2CPts" "N_Diff2Pts" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_InitRealMatrix" "N_KnotVectorGetKnots" "N_PtToCPt" "N_RealMatrixForBack" "N_RealMatrixInverse" "N_RealMatrixLuDecompose" "N_RealMatrixMultiply" "N_RealMatrixMultiplyCPtArray" "N_RealMatrixMultiplyPtArray" "N_RealMatrixTranspose" "N_ScaleCPt" "N_ScalePt" "N_SetRealMatrix" "N_VectorBlendCPt" "N_VectorBlendPt"))
 (:c "N_FitCrvApproxClosedConditions" :lisp "fit-crv-approx-closed-conditions" :file "NL_CrvFit.cpp" :line 14621 :lines 86 :wave 4 :status :done :root t
   :callees ("N_AllocKnotVectorAndArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_CrvDerivs" "N_CrvInitArrays" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCrvCalcKnotVector" "N_FitCrvHighDerivs" "N_FitCrvInterpGivenParams" "N_InitNurbs" "N_VectorCopy"))
 (:c "N_FitCrvHighDerivs" :lisp "fit-crv-high-derivs" :file "NL_CrvFit.cpp" :line 14829 :lines 101 :wave 4 :status :done
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisDerivs" "N_BasisEval" "N_CopyPt" "N_CrvGetCPtsKnotVectorAndKnots" "N_CrvSizeArrays" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCalcKnotsHighEndDerivs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_SetRealMatrix"))
 (:c "N_FitCalcKnotsDerivs" :lisp "fit-calc-knots-derivs" :file "NL_CrvFit.cpp" :line 15978 :lines 42 :wave 4 :status :done
   :callees ("N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_FitCalcKnotsHighEndDerivs" :lisp "fit-calc-knots-high-end-derivs" :file "NL_CrvFit.cpp" :line 16227 :lines 25 :wave 4 :status :done
   :callees ("N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_FitCrvPtsNormals" :lisp "fit-crv-pts-normals" :file "NL_CrvFit.cpp" :line 18067 :lines 44 :wave 4 :status :done :root t :fenced t
   :callees ("N_AllocPt1dArray" "N_CrvClosestPt" "N_CrvDerivs" "N_CrvGetParamBounds" "N_CrvInitArrays" "N_EndNurbs" "N_FitCrvFirstDeriv" "N_FitCrvInterp" "N_InitNurbs" "N_VectorCross"))
 (:c "ST_CalcCPtsCircle" :lisp "crv-fit--st-calc-c-pts-circle" :file "NL_CrvFit.cpp" :line 19297 :lines 16 :wave 4 :status :done :static t
   :callees ("N_DistPtPt" "N_Weight"))
 (:c "N_CrvClosestPt" :lisp "crv-closest-pt" :file "NL_CrvGeom.cpp" :line 64 :lines 100 :wave 2 :status :done
   :callees ("N_CopyPt" "N_CrvDerivs" "N_CrvGetKnots" "N_CrvIsClosed" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude"))
 (:c "N_CrvArcLength" :lisp "crv-arc-length" :file "NL_CrvGeom.cpp" :line 659 :lines 130 :wave 2 :status :done :root t
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_BezCalcLength" "N_Combine2CPts" "N_CopyCPt" "N_CrvExtractCrvSeg" "N_CrvGetBBoxMaxDiagDist" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvGetType" "N_CrvInitArrays" "N_DistPtPt" "N_EndNurbs" "N_InitNurbs" "N_VectorSum"))
 (:c "N_ApproxCrvWithPolyline" :lisp "approx-crv-with-polyline" :file "NL_CrvGeom.cpp" :line 1531 :lines 252 :wave 2 :status :done
   :callees ("N_AllocCPt2dArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BezCalcLength" "N_CPtToPtEuclid" "N_Combine2CPts" "N_CopyCPt" "N_CopyPt" "N_CrvDecomposeBez" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetDegree" "N_DistPtPt" "N_DistSqPtPt" "N_EndNurbs" "N_InitNurbs" "N_Realloc1dPtArray" "N_Realloc1dRealArray" "N_Realloc2dCPtArray" "N_Realloc2dRealArray" "N_VectorDir" "N_VectorDot" "N_VectorNormalize" "N_VectorPtAlongVector"))
 (:c "N_CrvPercentageAlongPt" :lisp "crv-percentage-along-pt" :file "NL_CrvGeom.cpp" :line 1945 :lines 204 :wave 2 :status :done :root t
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_BezCalcLength" "N_Combine2CPts" "N_CopyCPt" "N_CrvArcLength" "N_CrvClosestPt" "N_CrvEval" "N_CrvExtractCrvSeg" "N_CrvGetBBoxMaxDiagDist" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetDegree" "N_CrvGetKnotVector" "N_CrvGetKnots" "N_CrvGetType" "N_CrvInitArrays" "N_CrvReverse" "N_EndNurbs" "N_InitNurbs" "N_TranslateSum2Pts" "ST_bezlen"))
 (:c "N_CrvProjectPts" :lisp "crv-project-pts" :file "NL_CrvGeom.cpp" :line 2315 :lines 44 :wave 2 :status :done
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_ApproxCrvWithPolyline" "N_EndNurbs" "N_InitNurbs" "N_ProjectPtLineParam" "N_VectorCopy"))
 (:c "ST_bezlen" :lisp "crv-geom--st-bezlen" :file "NL_CrvGeom.cpp" :line 2390 :lines 56 :wave 2 :status :done :static t
   :callees ("N_AllocCPt1dArray" "N_BezSplit" "N_DistCPolygon" "N_DistCptCpt" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_CrvEvalCurvature" :lisp "crv-eval-curvature" :file "NL_CrvNurb.cpp" :line 663 :lines 40 :wave 0 :status :done :root t
   :callees ("N_CopyPt" "N_CrvDerivs" "N_CrvEvalFrenetFrame" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_VectorCross" "N_VectorMagnitude" "N_VectorNormalize" "N_VectorPtAlongVector"))
 (:c "N_CrvEval" :lisp "crv-eval" :file "NL_CrvNurb.cpp" :line 773 :lines 41 :wave 0 :status :done :root t
   :callees ("N_BasisEval" "N_CPtToPtEuclid" "N_KnotVectorIsParamOutOfBounds"))
 (:c "N_CrvGetBBox" :lisp "crv-get-b-box" :file "NL_CrvNurb.cpp" :line 873 :lines 16 :wave 0 :status :done
   :callees ("N_CrvGetArraySizes" "N_CrvGetEPolygon" "N_EndNurbs" "N_InitNurbs" "N_PolygonGetBBox"))
 (:c "N_CrvDerivs" :lisp "crv-derivs" :file "NL_CrvNurb.cpp" :line 955 :lines 177 :wave 0 :status :done
   :callees ("N_BasisFindSpan" "N_CPtToPtEuclid" "N_ErrSet" "N_FloatOpIsBad" "N_IsCrvRat" "N_SwapIntegers"))
 (:c "N_CrvEvalTangent" :lisp "crv-eval-tangent" :file "NL_CrvNurb.cpp" :line 1247 :lines 18 :wave 0 :status :done :root t
   :callees ("N_CopyPt" "N_CrvDerivs" "N_EndNurbs" "N_InitNurbs" "N_VectorNormalize"))
 (:c "N_CrvEvalFrenetFrame" :lisp "crv-eval-frenet-frame" :file "NL_CrvNurb.cpp" :line 1326 :lines 28 :wave 0 :status :done
   :callees ("N_CopyPt" "N_CrvDerivs" "N_CrvGetDegree" "N_EndNurbs" "N_InitNurbs" "N_VectorCross" "N_VectorCrossRef" "N_VectorNormalize" "N_VectorNormalizeRef"))
 (:c "N_CrvTransform" :lisp "crv-transform" :file "NL_CrvNurb.cpp" :line 1416 :lines 12 :wave 0 :status :done
   :callees ("N_Crv2dTo3d" "N_CrvGetCPts" "N_CrvIs3d" "N_TransformCPt"))
 (:c "N_CrvGetMinMaxWeightsAndPts" :lisp "crv-get-min-max-weights-and-pts" :file "NL_CrvNurb.cpp" :line 1480 :lines 44 :wave 0 :status :done
   :callees ("N_CPtGetW" "N_CPtToPtEuclid" "N_CrvGetCPts" "N_Diff2Pts" "N_IsCrvRat" "N_ScalePt" "N_Sum2Pts" "N_VectorMagnitude"))
 (:c "N_CrvLineFromPtAndVector" :lisp "crv-line-from-pt-and-vector" :file "NL_CrvNurb.cpp" :line 1883 :lines 21 :wave 0 :status :done :root t
   :callees ("N_CrvGetCPtsAndKnots" "N_CrvSizeArrays" "N_PtToCPt" "N_VectorSum"))
 (:c "N_CrvTranslate" :lisp "crv-translate" :file "NL_CrvNurb.cpp" :line 2083 :lines 15 :wave 0 :status :done
   :callees ("N_Crv2dTo3d" "N_CrvGetCPts" "N_CrvIs3d" "N_EndNurbs" "N_InitNurbs" "N_TranslateCPt"))
 (:c "N_CrvReverse" :lisp "crv-reverse" :file "NL_CrvNurb.cpp" :line 2430 :lines 36 :wave 0 :status :done :root t
   :callees ("N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvSizeArrays" "N_SwapCPts"))
 (:c "N_CrvGetBBoxMaxDiagDist" :lisp "crv-get-b-box-max-diag-dist" :file "NL_CrvNurb.cpp" :line 2521 :lines 14 :wave 0 :status :done :root t
   :callees ("N_BBoxGetDiagonal" "N_CrvGetBBox" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_CrvReparam" :lisp "crv-reparam" :file "NL_CrvNurb.cpp" :line 3454 :lines 6 :wave 0 :status :done :root t
   :callees ("N_CreateInterval" "N_CrvReparamToInterval"))
 (:c "N_CrvGetType" :lisp "crv-get-type" :file "NL_CrvNurb.cpp" :line 3531 :lines 229 :wave 0 :status :done
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_CPtToPtEuclid" "N_CopyPt" "N_CreatePlanePtNormal" "N_CrvEval" "N_CrvEvalCurvature" "N_CrvGetCPtsDegreeAndKnots" "N_Diff2Pts" "N_DistPtPt" "N_EndNurbs" "N_ErrClear" "N_FloatOpIsBad" "N_InitNurbs" "N_ProjectPtPlane" "N_ScalePt" "N_Sum2Pts" "N_VectorCross" "N_VectorDir" "N_VectorDot" "N_VectorMagnitude" "N_VectorNormalize" "N_VectorPtAlongVector" "N_VectorsAngle"))
 (:c "N_CrvGetCurvatureDeriv" :lisp "crv-get-curvature-deriv" :file "NL_CrvNurb.cpp" :line 3943 :lines 39 :wave 0 :status :done
   :callees ("N_CrvDerivs" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_VectorCross" "N_VectorDot" "N_VectorMagnitude" "N_VectorNormalizeRef"))
 (:c "N_CrvOffsetPtSampling" :lisp "crv-offset-pt-sampling" :file "NL_CrvOffset.cpp" :line 1217 :lines 213 :wave 5 :status :done
   :callees ("N_AllocInt1dArray" "N_AllocKnotVectorAndArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_CreateCircArc" "N_CreateInterval" "N_Crv3dTo2d" "N_CrvCompress" "N_CrvDecomposeBez" "N_CrvDerivs" "N_CrvEvalTangent" "N_CrvGetArraySizes" "N_CrvGetDegree" "N_CrvGetKnots" "N_CrvGetType" "N_CrvInitArrays" "N_CrvIs3d" "N_CrvLineFromPtAndVector" "N_CrvOffsetGetMax2ndDeriv" "N_CrvRemoveKnotsDerivConstraints" "N_CrvReparamToInterval" "N_DistPolygon" "N_EndNurbs" "N_FitCalcCrvParamValues" "N_FitCalcKnotVectorEndDerivs" "N_FitCrvCalcKnotVector" "N_FitCrvInterpGivenParams" "N_FitCrvKnotsAndDerivs" "N_FloatOpIsBad" "N_InitNurbs" "N_VectorBlendPt" "N_VectorCopy" "N_VectorCross" "N_VectorNormalizeRef" "N_VectorPtAlongVector" "N_VectorScale"))
 (:c "N_CrvInsertKnot" :lisp "crv-insert-knot" :file "NL_CrvTool.cpp" :line 65 :lines 79 :wave 1 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocCrvArrays" "N_BasisFindSpanAndMult" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvSizeArrays" "N_EndNurbs" "N_FreeCrv" "N_InitNurbs" "N_KnotVectorIsEndParam"))
 (:c "N_CrvDecomposeBez" :lisp "crv-decompose-bez" :file "NL_CrvTool.cpp" :line 811 :lines 109 :wave 1 :status :done :root t
   :callees ("N_Alloc1dArrayCrvs" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPts" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_CrvDecomposeContinuity" :lisp "crv-decompose-continuity" :file "NL_CrvTool.cpp" :line 1031 :lines 74 :wave 1 :status :done :root t
   :callees ("N_AllocArrayCrvPtrs" "N_AllocCrvAndArrays" "N_AllocInt1dArray" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_CrvRefine" :lisp "crv-refine" :file "NL_CrvTool.cpp" :line 1369 :lines 100 :wave 1 :status :done
   :callees ("N_AllocCrvArrays" "N_BasisFindSpan" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvSizeArrays" "N_EndNurbs" "N_FreeCrv" "N_InitNurbs" "N_KnotVectorGetKnots" "N_KnotVectorIsEndParam"))
 (:c "N_CrvExtractCrvSeg" :lisp "crv-extract-crv-seg" :file "NL_CrvTool.cpp" :line 1561 :lines 86 :wave 1 :status :done
   :callees ("N_AllocCrvArrays" "N_BasisFindSpanAndMult" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvSizeArrays" "N_FreeCrv"))
 (:c "N_CrvRemoveKnot" :lisp "crv-remove-knot" :file "NL_CrvTool.cpp" :line 1752 :lines 175 :wave 1 :status :done
   :callees ("N_AllocCPt1dArray" "N_BasisFindGlobalMax" "N_BasisFindSpanAndMult" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvGetMinMaxWeightsAndPts" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_DistCptCptHomo" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat" "N_KnotVectorIsEndParam"))
 (:c "N_CrvRemoveKnotMaxErr" :lisp "crv-remove-knot-max-err" :file "NL_CrvTool.cpp" :line 2174 :lines 47 :wave 1 :status :done
   :callees ("N_AllocCPt1dArray" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsDegreeAndKnots" "N_DistCptCptHomo" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_CrvRemoveKnots" :lisp "crv-remove-knots" :file "NL_CrvTool.cpp" :line 2295 :lines 279 :wave 1 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocInt1dArray" "N_AllocReal1dArray" "N_BasisFindAllSpanMaxima" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_CrvCompress" "N_CrvCopy" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetCPtsKnotVectorAndKnots" "N_CrvGetMinMaxWeightsAndPts" "N_CrvRemoveKnotMaxErr" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat"))
 (:c "N_CrvElevateDegree" :lisp "crv-elevate-degree" :file "NL_CrvTool.cpp" :line 3285 :lines 191 :wave 1 :status :done :root t
   :callees ("N_AllocCPt1dArray" "N_AllocCrvArrays" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_BezElevateDegree" "N_BezGetDegreeElevationMatrix" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvSizeArrays" "N_EndNurbs" "N_FreeCrv" "N_InitNurbs" "N_InitRealMatrix" "N_KnotVectorIsValid"))
 (:c "N_CrvReduceDegreeOnce" :lisp "crv-reduce-degree-once" :file "NL_CrvTool.cpp" :line 3636 :lines 352 :wave 1 :status :done :root t
   :callees ("N_AllocCPt1dArray" "N_AllocCrvArrays" "N_AllocReal1dArray" "N_BasisFindAllSpanMaxima" "N_BasisGetSpanCount" "N_BezDegreeReduceCoefs" "N_BezReduceDegree" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvGetMinMaxWeightsAndPts" "N_CrvInitArrays" "N_CrvSizeArrays" "N_DistCptCptHomo" "N_EndNurbs" "N_FreeCrv" "N_InitNurbs" "N_IsCrvRat" "N_KnotVectorFromRealArray"))
 (:c "N_CrvsMakeCompatibleKnotTol" :lisp "crvs-make-compatible-knot-tol" :file "NL_CrvTool.cpp" :line 4289 :lines 81 :wave 1 :status :done
   :callees ("N_Alloc1dArrayKnotVectPtrs" "N_CrvElevateDegree" "N_CrvGetDegree" "N_CrvGetKnotVector" "N_CrvGetParamBounds" "N_CrvRefine" "N_CrvReparamToInterval" "N_CrvsMakeRatCompatible" "N_EndNurbs" "N_GetCompatibleKnotArray" "N_GetCompatibleKnotVectorToTol" "N_InitNurbs" "N_KnotVectorGetKnots"))
 (:c "N_CrvsMakeCompatible" :lisp "crvs-make-compatible" :file "NL_CrvTool.cpp" :line 4468 :lines 5 :wave 1 :status :done :root t
   :callees ("N_CrvsMakeCompatibleKnotTol"))
 (:c "N_AllocCPolygon" :lisp "alloc-c-polygon" :file "NL_FrameBasic.cpp" :line 42 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocCPolygonAndArray" :lisp "alloc-c-polygon-and-array" :file "NL_FrameBasic.cpp" :line 112 :lines 13 :wave 0 :status :vanish
   :callees ("N_AllocCPolygon" "N_AllocCPt1dArray" "N_CPolygonFromCPts"))
 (:c "N_FreeCPolygon" :lisp "free-c-polygon" :file "NL_FrameBasic.cpp" :line 164 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_CPolygonFromCPts" :lisp "c-polygon-from-c-pts" :file "NL_FrameBasic.cpp" :line 276 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_PolygonGetBBox" :lisp "polygon-get-b-box" :file "NL_FrameBasic.cpp" :line 795 :lines 28 :wave 0 :status :done
   :callees ("N_BBoxDefine" "N_EPolygonGetPts" "N_PtToXYZ"))
 (:c "N_EPolygonFromPts" :lisp "e-polygon-from-pts" :file "NL_FrameBasic.cpp" :line 964 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_EPolygonGetPts" :lisp "e-polygon-get-pts" :file "NL_FrameBasic.cpp" :line 1002 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_AllocCNet" :lisp "alloc-c-net" :file "NL_FrameBasic.cpp" :line 1734 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocCNetAndArrays" :lisp "alloc-c-net-and-arrays" :file "NL_FrameBasic.cpp" :line 1803 :lines 13 :wave 0 :status :vanish
   :callees ("N_AllocCNet" "N_AllocCPt2dArray" "N_CNetFromCPts"))
 (:c "N_FreeCNet" :lisp "free-c-net" :file "NL_FrameBasic.cpp" :line 1859 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_CNetFromCPts" :lisp "c-net-from-c-pts" :file "NL_FrameBasic.cpp" :line 1936 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_ENetFromPts" :lisp "e-net-from-pts" :file "NL_FrameBasic.cpp" :line 2082 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_ENetGetBBox" :lisp "e-net-get-b-box" :file "NL_FrameBasic.cpp" :line 2117 :lines 31 :wave 0 :status :done
   :callees ("N_BBoxDefine" "N_ENetGetPts" "N_PtToXYZ"))
 (:c "N_ENetGetPts" :lisp "e-net-get-pts" :file "NL_FrameBasic.cpp" :line 2196 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_ENetIsClosed" :lisp "e-net-is-closed" :file "NL_FrameBasic.cpp" :line 2234 :lines 36 :wave 0 :status :done
   :callees ("N_DistPtPt" "N_ENetGetPts"))
 (:c "N_CrvFuncEvalInvertPt" :lisp "crv-func-eval-invert-pt" :file "NL_FuncsAdv.cpp" :line 264 :lines 75 :wave 2 :status :done
   :callees ("N_CFuncDerivs" "N_CFuncEval" "N_CFuncGetData" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs"))
 (:c "N_MapKnotsBetweenCrvFuncAndKnotVector" :lisp "map-knots-between-crv-func-and-knot-vector" :file "NL_FuncsAdv.cpp" :line 456 :lines 159 :wave 2 :status :done
   :callees ("N_BasisFindSpan" "N_BasisGetLongestAndShortestSpans" "N_BasisGetSpanCount" "N_CFuncEval" "N_CFuncGetData" "N_CrvFuncCntrlValKnotVector" "N_CrvFuncEvalInvertPt" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SetKnotIndex"))
 (:c "N_CrvFuncRefine" :lisp "crv-func-refine" :file "NL_FuncsAdv.cpp" :line 807 :lines 88 :wave 2 :status :done
   :callees ("N_AllocCFuncArrays" "N_BasisFindSpan" "N_CFuncGetData" "N_CFuncGetKnotVector" "N_CFuncSizeArrays" "N_CrvFuncCntrlValKnots" "N_EndNurbs" "N_FreeCrvFunc" "N_InitNurbs" "N_KnotVectorGetKnots" "N_KnotVectorIsEndParam"))
 (:c "N_CrvFuncDegreeElevate" :lisp "crv-func-degree-elevate" :file "NL_FuncsAdv.cpp" :line 1002 :lines 180 :wave 2 :status :done :fenced t
   :callees ("N_AllocCFuncArrays" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_BezFuncDegreeElevate" "N_BezGetDegreeElevationMatrix" "N_CFuncGetData" "N_CFuncGetKnotVector" "N_CFuncSizeArrays" "N_CrvFuncCntrlValKnots" "N_EndNurbs" "N_FreeCrvFunc" "N_InitNurbs" "N_InitRealMatrix" "N_KnotVectorIsValid"))
 (:c "N_CrvFuncMultiplyCrvFunc" :lisp "crv-func-multiply-crv-func" :file "NL_FuncsAdv.cpp" :line 1705 :lines 328 :wave 2 :status :done
   :callees ("N_AllocKnotVectorAndArray" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_BezFuncMultiplyBezCrvMatrix" "N_BezFuncMultiplyBezFunc" "N_CFuncGetArraySizes" "N_CFuncGetData" "N_CFuncGetDegree" "N_CFuncGetKnotVector" "N_CFuncInitArrays" "N_CFuncSetSizeIndices" "N_CFuncSizeArrays" "N_CrvFuncCntrlValKnots" "N_CrvFuncCompact" "N_CrvFuncRefine" "N_EndNurbs" "N_InitNurbs" "N_InitRealMatrix" "N_KnotVectorGetKnots" "N_MakeKnotsCompatible"))
 (:c "N_CrvFuncMultiplyCrv" :lisp "crv-func-multiply-crv" :file "NL_FuncsAdv.cpp" :line 2229 :lines 330 :wave 2 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocKnotVectorAndArray" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_BezFuncMultiplyBezCrv" "N_BezFuncMultiplyBezCrvMatrix" "N_CFuncGetArraySizes" "N_CFuncGetData" "N_CFuncGetDegree" "N_CFuncGetKnotVector" "N_CFuncInitArrays" "N_Combine2CPts" "N_CopyCPt" "N_CrvCompress" "N_CrvFuncCntrlValKnots" "N_CrvFuncRefine" "N_CrvGetArraySizes" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetDegree" "N_CrvGetKnotVector" "N_CrvInitArrays" "N_CrvRefine" "N_CrvSetSizeIndices" "N_CrvSizeArrays" "N_EndNurbs" "N_InitNurbs" "N_InitRealMatrix" "N_KnotVectorGetKnots" "N_MakeKnotsCompatible" "N_ScaleCPtXYZ"))
 (:c "N_CrvFuncSumDiffCrvFunc" :lisp "crv-func-sum-diff-crv-func" :file "NL_FuncsAdv.cpp" :line 2767 :lines 102 :wave 2 :status :done :fenced t
   :callees ("N_AllocKnotVectorAndArray" "N_BasisGetSpanCount" "N_CFuncGetData" "N_CFuncGetKnotVector" "N_CFuncInitArrays" "N_CFuncSizeArrays" "N_CrvFuncCntrlValKnots" "N_CrvFuncCopy" "N_CrvFuncDegreeElevate" "N_CrvFuncRefine" "N_EndNurbs" "N_GetCompatibleKnotArrayMult" "N_InitNurbs" "N_KnotVectorGetKnots"))
 (:c "N_CrvFuncDeriv" :lisp "crv-func-deriv" :file "NL_FuncsAdv.cpp" :line 2975 :lines 49 :wave 2 :status :done
   :callees ("N_AllocReal1dArray" "N_CFuncGetData" "N_CFuncSizeArrays" "N_CrvFuncCntrlValKnots" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs"))
 (:c "N_CrvFuncMultiplyConstant" :lisp "crv-func-multiply-constant" :file "NL_FuncsAdv.cpp" :line 3975 :lines 8 :wave 2 :status :done :fenced t
   :callees ("N_CrvFuncCntrlVal"))
 (:c "N_CrvIsFuncSized" :lisp "crv-is-func-sized" :file "NL_FuncsBasic.cpp" :line 47 :lines 12 :wave 0 :status :done
   :callees ("N_CFuncGetArraySizes" "N_ErrSet"))
 (:c "N_AllocCValue" :lisp "alloc-c-value" :file "NL_FuncsBasic.cpp" :line 93 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocCValueAndArray" :lisp "alloc-c-value-and-array" :file "NL_FuncsBasic.cpp" :line 163 :lines 13 :wave 0 :status :vanish
   :callees ("N_AllocCValue" "N_AllocReal1dArray" "N_CValueFromArray"))
 (:c "N_FreeCValue" :lisp "free-c-value" :file "NL_FuncsBasic.cpp" :line 219 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_CValueFromArray" :lisp "c-value-from-array" :file "NL_FuncsBasic.cpp" :line 296 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_AllocCrvFuncArray" :lisp "alloc-crv-func-array" :file "NL_FuncsBasic.cpp" :line 339 :lines 35 :wave 0 :status :vanish
   :callees ("N_AllocCrvFunc" "N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocCrvFunc" :lisp "alloc-crv-func" :file "NL_FuncsBasic.cpp" :line 428 :lines 17 :wave 0 :status :vanish
   :callees ("N_AllocCValueAndArray" "N_AllocCrvFuncStack" "N_AllocKnotVectorAndArray" "N_CFuncFromKnotVector"))
 (:c "N_AllocCrvFuncStack" :lisp "alloc-crv-func-stack" :file "NL_FuncsBasic.cpp" :line 489 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_FreeCrvFunc" :lisp "free-crv-func" :file "NL_FuncsBasic.cpp" :line 556 :lines 13 :wave 0 :status :vanish
   :callees ("N_CFuncGetArrayAndKnotVector" "N_CrvFuncCntrlValKnots" "N_FreeCValue" "N_FreeKnotVector" "N_FreeReal1dArray"))
 (:c "N_CFuncFromKnotVector" :lisp "c-func-from-knot-vector" :file "NL_FuncsBasic.cpp" :line 613 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_CFuncFromKnots" :lisp "c-func-from-knots" :file "NL_FuncsBasic.cpp" :line 661 :lines 15 :wave 0 :status :done
   :callees ("N_AllocCValue" "N_AllocKnotVector" "N_CFuncFromKnotVector" "N_CValueFromArray" "N_KnotVectorFromRealArray"))
 (:c "N_CrvAreFuncArraysNULL" :lisp "crv-are-func-arrays-null" :file "NL_FuncsBasic.cpp" :line 717 :lines 15 :wave 0 :status :vanish
   :callees ("N_CFuncGetArrayAndKnotVector"))
 (:c "N_AllocCFuncArrays" :lisp "alloc-c-func-arrays" :file "NL_FuncsBasic.cpp" :line 780 :lines 13 :wave 0 :status :vanish
   :callees ("N_AllocCValueAndArray" "N_AllocKnotVectorAndArray" "N_CFuncFromKnotVector"))
 (:c "N_CFuncInitArrays" :lisp "c-func-init-arrays" :file "NL_FuncsBasic.cpp" :line 837 :lines 6 :wave 0 :status :vanish
   :callees ())
 (:c "N_CFuncSetSizeIndices" :lisp "c-func-set-size-indices" :file "NL_FuncsBasic.cpp" :line 874 :lines 6 :wave 0 :status :vanish
   :callees ())
 (:c "N_CFuncSizeArrays" :lisp "c-func-size-arrays" :file "NL_FuncsBasic.cpp" :line 922 :lines 18 :wave 0 :status :vanish
   :callees ("N_AllocCFuncArrays" "N_CFuncSetSizeIndices" "N_CrvAreFuncArraysNULL" "N_CrvIsFuncSized"))
 (:c "N_CrvFuncCompact" :lisp "crv-func-compact" :file "NL_FuncsBasic.cpp" :line 985 :lines 23 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_CFuncGetArraySizes" "N_CrvFuncCntrlValKnots" "N_CrvFuncSetPtrs" "N_FreeReal1dArray"))
 (:c "N_CrvFuncCopy" :lisp "crv-func-copy" :file "NL_FuncsBasic.cpp" :line 1066 :lines 18 :wave 0 :status :done
   :callees ("N_CFuncGetData" "N_CFuncSizeArrays" "N_CrvFuncCntrlValKnots"))
 (:c "N_CrvFuncCntrlVal" :lisp "crv-func-cntrl-val" :file "NL_FuncsBasic.cpp" :line 1273 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_CFuncGetArrayAndKnotVector" :lisp "c-func-get-array-and-knot-vector" :file "NL_FuncsBasic.cpp" :line 1311 :lines 6 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvFuncCntrlValKnotVector" :lisp "crv-func-cntrl-val-knot-vector" :file "NL_FuncsBasic.cpp" :line 1348 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_CFuncGetArraySizes" :lisp "c-func-get-array-sizes" :file "NL_FuncsBasic.cpp" :line 1382 :lines 5 :wave 0 :status :vanish
   :callees ())
 (:c "N_CFuncGetKnots" :lisp "c-func-get-knots" :file "NL_FuncsBasic.cpp" :line 1417 :lines 5 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_CFuncGetKnotVector" :lisp "c-func-get-knot-vector" :file "NL_FuncsBasic.cpp" :line 1450 :lines 4 :wave 0 :status :vanish
   :callees ())
 (:c "N_CFuncGetDegree" :lisp "c-func-get-degree" :file "NL_FuncsBasic.cpp" :line 1481 :lines 4 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvFuncCntrlValKnots" :lisp "crv-func-cntrl-val-knots" :file "NL_FuncsBasic.cpp" :line 1513 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_CFuncGetData" :lisp "c-func-get-data" :file "NL_FuncsBasic.cpp" :line 1552 :lines 8 :wave 0 :status :vanish
   :callees ())
 (:c "N_CrvFuncSetPtrs" :lisp "crv-func-set-ptrs" :file "NL_FuncsBasic.cpp" :line 1588 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_CFuncEval" :lisp "c-func-eval" :file "NL_FuncsBasic.cpp" :line 1630 :lines 31 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_BasisEval" "N_CFuncGetDegree" "N_CFuncGetKnotVector" "N_CrvFuncCntrlValKnots" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorIsParamOutOfBounds"))
 (:c "N_CFuncDerivs" :lisp "c-func-derivs" :file "NL_FuncsBasic.cpp" :line 1736 :lines 34 :wave 0 :status :done
   :callees ("N_AllocReal2dArray" "N_BasisDerivs" "N_CFuncGetDegree" "N_CFuncGetKnotVector" "N_CrvFuncCntrlValKnots" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorIsParamOutOfBounds"))
 (:c "N_SrfIsFuncSized" :lisp "srf-is-func-sized" :file "NL_FuncsBasic.cpp" :line 1833 :lines 12 :wave 0 :status :done
   :callees ("N_ErrSet" "N_SFuncGetArraySizes"))
 (:c "N_AllocSValue" :lisp "alloc-s-value" :file "NL_FuncsBasic.cpp" :line 1883 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocSValueAndArray" :lisp "alloc-s-value-and-array" :file "NL_FuncsBasic.cpp" :line 1953 :lines 13 :wave 0 :status :vanish
   :callees ("N_AllocReal2dArray" "N_AllocSValue" "N_SValueFromArray"))
 (:c "N_SValueFromArray" :lisp "s-value-from-array" :file "NL_FuncsBasic.cpp" :line 2086 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_SFuncFromKnotVectors" :lisp "s-func-from-knot-vectors" :file "NL_FuncsBasic.cpp" :line 2460 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_SrfAreFuncArraysNULL" :lisp "srf-are-func-arrays-null" :file "NL_FuncsBasic.cpp" :line 2572 :lines 15 :wave 0 :status :vanish
   :callees ("N_SFuncGetArrayAndKnotVectors"))
 (:c "N_AllocSFuncArrays" :lisp "alloc-s-func-arrays" :file "NL_FuncsBasic.cpp" :line 2631 :lines 16 :wave 0 :status :vanish
   :callees ("N_AllocKnotVectorAndArray" "N_AllocSValueAndArray" "N_SFuncFromKnotVectors"))
 (:c "N_SFuncSetSizeIndices" :lisp "s-func-set-size-indices" :file "NL_FuncsBasic.cpp" :line 2730 :lines 9 :wave 0 :status :vanish
   :callees ())
 (:c "N_SFuncSizeArrays" :lisp "s-func-size-arrays" :file "NL_FuncsBasic.cpp" :line 2779 :lines 18 :wave 0 :status :vanish
   :callees ("N_AllocSFuncArrays" "N_SFuncSetSizeIndices" "N_SrfAreFuncArraysNULL" "N_SrfIsFuncSized"))
 (:c "N_SrfFuncCntrlVal" :lisp "srf-func-cntrl-val" :file "NL_FuncsBasic.cpp" :line 3186 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_SFuncGetKnotVectors" :lisp "s-func-get-knot-vectors" :file "NL_FuncsBasic.cpp" :line 3343 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_SFuncGetArraySizes" :lisp "s-func-get-array-sizes" :file "NL_FuncsBasic.cpp" :line 3408 :lines 7 :wave 0 :status :vanish
   :callees ())
 (:c "N_SFuncGetArrayAndKnotVectors" :lisp "s-func-get-array-and-knot-vectors" :file "NL_FuncsBasic.cpp" :line 3446 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_VectorDot" :lisp "vector-dot" :file "NL_Geom.cpp" :line 51 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_VectorDotRef" :lisp "vector-dot-ref" :file "NL_Geom.cpp" :line 84 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_VectorCross" :lisp "vector-cross" :file "NL_Geom.cpp" :line 116 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorCrossRef" :lisp "vector-cross-ref" :file "NL_Geom.cpp" :line 151 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorScale" :lisp "vector-scale" :file "NL_Geom.cpp" :line 188 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorReverse" :lisp "vector-reverse" :file "NL_Geom.cpp" :line 222 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorReverseInPlace" :lisp "vector-reverse-in-place" :file "NL_Geom.cpp" :line 256 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorScaleRef" :lisp "vector-scale-ref" :file "NL_Geom.cpp" :line 375 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorSum" :lisp "vector-sum" :file "NL_Geom.cpp" :line 410 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorCombine" :lisp "vector-combine" :file "NL_Geom.cpp" :line 448 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorDiff" :lisp "vector-diff" :file "NL_Geom.cpp" :line 483 :lines 9 :wave 0 :status :done
   :callees ())
 (:c "N_VectorDiffCPts" :lisp "vector-diff-c-pts" :file "NL_Geom.cpp" :line 524 :lines 7 :wave 0 :status :port
   :callees ("N_CPtToPtEuclid" "N_VectorDiff"))
 (:c "N_VectorDir" :lisp "vector-dir" :file "NL_Geom.cpp" :line 566 :lines 9 :wave 0 :status :done
   :callees ())
 (:c "N_VectorCreate" :lisp "vector-create" :file "NL_Geom.cpp" :line 604 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorCopy" :lisp "vector-copy" :file "NL_Geom.cpp" :line 639 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_VectorMagnitude" :lisp "vector-magnitude" :file "NL_Geom.cpp" :line 714 :lines 6 :wave 0 :status :done
   :callees ("N_VectorDot"))
 (:c "N_VectorMagnitudeRef" :lisp "vector-magnitude-ref" :file "NL_Geom.cpp" :line 752 :lines 6 :wave 0 :status :done
   :callees ("N_VectorDotRef"))
 (:c "N_VectorNormalize" :lisp "vector-normalize" :file "NL_Geom.cpp" :line 791 :lines 20 :wave 0 :status :done
   :callees ("N_ErrSet" "N_FloatOpIsBad" "N_VectorMagnitude" "N_VectorScale"))
 (:c "N_VectorNormalizeRef" :lisp "vector-normalize-ref" :file "NL_Geom.cpp" :line 843 :lines 14 :wave 0 :status :done
   :callees ("N_ErrSet" "N_FloatOpIsBad" "N_VectorMagnitudeRef" "N_VectorScaleRef"))
 (:c "N_VectorPerpendicular" :lisp "vector-perpendicular" :file "NL_Geom.cpp" :line 897 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_VectorPtAlongVector" :lisp "vector-pt-along-vector" :file "NL_Geom.cpp" :line 940 :lines 10 :wave 0 :status :done
   :callees ())
 (:c "N_VectorsAngle" :lisp "vectors-angle" :file "NL_Geom.cpp" :line 1077 :lines 21 :wave 0 :status :done
   :callees ("N_ErrSet" "N_FloatOpIsBad" "N_VectorDot" "N_VectorMagnitude"))
 (:c "N_VectorDirectedAngle" :lisp "vector-directed-angle" :file "NL_Geom.cpp" :line 1137 :lines 27 :wave 0 :status :done
   :callees ("N_ErrSet" "N_FloatOpIsBad" "N_VectorDot" "N_VectorMagnitude" "N_VectorPerpendicular"))
 (:c "N_IsectLineLine" :lisp "isect-line-line" :file "NL_Geom.cpp" :line 1215 :lines 67 :wave 0 :status :done
   :callees ("N_Combine2Pts" "N_ErrSet" "N_FloatOpIsBad" "N_LineGetData" "N_VectorDir" "N_VectorDot" "N_VectorNormalize" "N_VectorPtAlongVector"))
 (:c "N_LineSegs2dAreIntersecting" :lisp "line-segs-2d-are-intersecting" :file "NL_Geom.cpp" :line 1449 :lines 46 :wave 0 :status :port
   :callees ())
 (:c "N_IsectLinePlane" :lisp "isect-line-plane" :file "NL_Geom.cpp" :line 1552 :lines 47 :wave 0 :status :done
   :callees ("N_ErrSet" "N_FloatOpIsBad" "N_LineGetData" "N_PlaneGetData" "N_VectorDir" "N_VectorDot" "N_VectorNormalize" "N_VectorPtAlongVector"))
 (:c "N_IsectPlanePlane" :lisp "isect-plane-plane" :file "NL_Geom.cpp" :line 1662 :lines 35 :wave 0 :status :done
   :callees ("N_PlaneGetData" "N_PlanePtNormalToImplicit" "N_PtFromXYZ"))
 (:c "N_ProjectPtLine" :lisp "project-pt-line" :file "NL_Geom.cpp" :line 1757 :lines 8 :wave 0 :status :done
   :callees ("N_CreatePlanePtNormal" "N_IsectLinePlane"))
 (:c "N_ProjectPtLineParam" :lisp "project-pt-line-param" :file "NL_Geom.cpp" :line 1809 :lines 55 :wave 0 :status :done
   :callees ("N_DistPtPt" "N_ErrSet" "N_FloatOpIsBad" "N_VectorCopy" "N_VectorDiff" "N_VectorDot" "N_VectorPtAlongVector" "N_VectorScale" "N_VectorSum"))
 (:c "N_ProjectPtPlane" :lisp "project-pt-plane" :file "NL_Geom.cpp" :line 1922 :lines 11 :wave 0 :status :done
   :callees ("N_CreateLineStartDirVector" "N_IsectLinePlane" "N_PlaneGetNormal"))
 (:c "N_DistPtPt" :lisp "dist-pt-pt" :file "NL_Geom.cpp" :line 2678 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_DistSqPtPt" :lisp "dist-sq-pt-pt" :file "NL_Geom.cpp" :line 2756 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_DistSqPtPt2d" :lisp "dist-sq-pt-pt-2d" :file "NL_Geom.cpp" :line 2790 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_PtsAreEqual" :lisp "pts-are-equal" :file "NL_Geom.cpp" :line 2823 :lines 11 :wave 0 :status :done
   :callees ())
 (:c "N_DistPtLineSeg" :lisp "dist-pt-line-seg" :file "NL_Geom.cpp" :line 2868 :lines 44 :wave 0 :status :done
   :callees ("N_DistPtPt" "N_ErrSet" "N_FloatOpIsBad" "N_VectorCross" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude" "N_VectorScale" "N_VectorSum"))
 (:c "N_DistPerpPtLineSeg" :lisp "dist-perp-pt-line-seg" :file "NL_Geom.cpp" :line 2963 :lines 37 :wave 0 :status :done
   :callees ("N_DistPtPt" "N_ErrSet" "N_FloatOpIsBad" "N_VectorCross" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude" "N_VectorScale" "N_VectorSum"))
 (:c "N_DistPolygon" :lisp "dist-polygon" :file "NL_Geom.cpp" :line 3050 :lines 11 :wave 0 :status :done
   :callees ("N_DistPtPt"))
 (:c "N_DistCPolygon" :lisp "dist-c-polygon" :file "NL_Geom.cpp" :line 3094 :lines 11 :wave 0 :status :done
   :callees ("N_DistCptCpt"))
 (:c "N_DistCptCpt" :lisp "dist-cpt-cpt" :file "NL_Geom.cpp" :line 3138 :lines 7 :wave 0 :status :done
   :callees ("N_CPtToPtEuclid" "N_DistPtPt"))
 (:c "N_DistCptCptHomo" :lisp "dist-cpt-cpt-homo" :file "NL_Geom.cpp" :line 3181 :lines 10 :wave 0 :status :done
   :callees ())
 (:c "N_DistPtInfLine" :lisp "dist-pt-inf-line" :file "NL_Geom.cpp" :line 3227 :lines 12 :wave 0 :status :done
   :callees ("N_DistPtPt" "N_LineSetBoundingFlag" "N_ProjectPtLine"))
 (:c "N_DistSignedPtPlane" :lisp "dist-signed-pt-plane" :file "NL_Geom.cpp" :line 3340 :lines 6 :wave 0 :status :done
   :callees ("N_PtToXYZ"))
 (:c "N_LineSetBoundingFlag" :lisp "line-set-bounding-flag" :file "NL_Geom.cpp" :line 3451 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_LineGetData" :lisp "line-get-data" :file "NL_Geom.cpp" :line 3523 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_Pts1dGetMaxExtent" :lisp "pts-1d-get-max-extent" :file "NL_Geom.cpp" :line 3560 :lines 8 :wave 0 :status :done
   :callees ("N_BBoxGetDiagonal" "N_EPolygonFromPts" "N_PolygonGetBBox"))
 (:c "N_Pts2dGetMaxExtent" :lisp "pts-2d-get-max-extent" :file "NL_Geom.cpp" :line 3603 :lines 8 :wave 0 :status :done
   :callees ("N_BBoxGetDiagonal" "N_ENetFromPts" "N_ENetGetBBox"))
 (:c "N_BBoxGetDiagonal" :lisp "b-box-get-diagonal" :file "NL_Geom.cpp" :line 3644 :lines 9 :wave 0 :status :done
   :callees ("N_DistPtPt" "N_GetBBoxData" "N_PtFromXYZ"))
 (:c "N_BBox2dCalcOverlap" :lisp "b-box-2d-calc-overlap" :file "NL_Geom.cpp" :line 3780 :lines 12 :wave 0 :status :port
   :callees ())
 (:c "N_CreateLinePtPt" :lisp "create-line-pt-pt" :file "NL_Geom.cpp" :line 3832 :lines 6 :wave 0 :status :done
   :callees ("N_CreateLineStartDirVector" "N_VectorDir"))
 (:c "N_CreateLineStartDirVector" :lisp "create-line-start-dir-vector" :file "NL_Geom.cpp" :line 3875 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_CreatePlanePtNormal" :lisp "create-plane-pt-normal" :file "NL_Geom.cpp" :line 3912 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_CreateInterval" :lisp "create-interval" :file "NL_Geom.cpp" :line 3951 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_IntervalGetData" :lisp "interval-get-data" :file "NL_Geom.cpp" :line 3987 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_CreateRectangle" :lisp "create-rectangle" :file "NL_Geom.cpp" :line 4025 :lines 7 :wave 0 :status :done
   :callees ())
 (:c "N_RectangleGetData" :lisp "rectangle-get-data" :file "NL_Geom.cpp" :line 4065 :lines 7 :wave 0 :status :done
   :callees ())
 (:c "N_GetBBoxData" :lisp "get-b-box-data" :file "NL_Geom.cpp" :line 4103 :lines 9 :wave 0 :status :done
   :callees ())
 (:c "N_BBoxDefine" :lisp "b-box-define" :file "NL_Geom.cpp" :line 4143 :lines 9 :wave 0 :status :done
   :callees ())
 (:c "N_PlaneGetData" :lisp "plane-get-data" :file "NL_Geom.cpp" :line 4183 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_PlaneGetNormal" :lisp "plane-get-normal" :file "NL_Geom.cpp" :line 4217 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_CreateTransformMatrixFromAxes" :lisp "create-transform-matrix-from-axes" :file "NL_Geom.cpp" :line 4387 :lines 97 :wave 0 :status :done
   :callees ("N_CheckMemRealMatrix" "N_EndNurbs" "N_GetRealMatrixPtr" "N_InitNurbs" "N_PtToXYZ" "N_RealMatrixMultiply" "N_SetRealMatrix"))
 (:c "N_CptReflect" :lisp "cpt-reflect" :file "NL_Geom.cpp" :line 5983 :lines 33 :wave 0 :status :done
   :callees ("N_CPtFromWxWyWz" "N_CPtToPtEuclid" "N_CPtToWxWyWz" "N_PlaneGetNormal" "N_ProjectPtPlane" "N_PtToXYZ" "N_Sum2Pts" "N_VectorDir" "N_VectorNormalize"))
 (:c "N_PlanePtNormalToImplicit" :lisp "plane-pt-normal-to-implicit" :file "NL_Geom.cpp" :line 6398 :lines 13 :wave 0 :status :done
   :callees ("N_PtToXYZ" "N_VectorDot" "N_VectorNormalizeRef"))
 (:c "N_DistPtArc" :lisp "dist-pt-arc" :file "NL_Geom.cpp" :line 6508 :lines 76 :wave 0 :status :done
   :callees ("N_DistPerpPtLineSeg" "N_DistPtPt" "N_VectorCreate" "N_VectorDiff" "N_VectorDirectedAngle"))
 (:c "N_DistMaxLineArc" :lisp "dist-max-line-arc" :file "NL_Geom.cpp" :line 6648 :lines 99 :wave 0 :status :done
   :callees ("N_CreateLinePtPt" "N_DistPerpPtLineSeg" "N_DistPtArc" "N_DistPtPt" "N_IsectLineLine" "N_ProjectPtLineParam" "N_VectorCopy"))
 (:c "N_LineImplicitToPtVector" :lisp "line-implicit-to-pt-vector" :file "NL_Geom.cpp" :line 6878 :lines 29 :wave 0 :status :done
   :callees ("N_FloatOpIsBad" "N_PtFromXYZ" "N_VectorNormalizeRef"))
 (:c "N_LineFit2dPts" :lisp "line-fit-2d-pts" :file "NL_Geom.cpp" :line 7030 :lines 148 :wave 0 :status :done
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_EndNurbs" "N_FloatOpIsBad" "N_GetRealMatrixPtr" "N_InitNurbs" "N_InitRealMatrix" "N_PtToXYZ" "N_RealMatrixInverse" "N_RealMatrixMultiply" "N_RealMatrixSingleValueDecompose" "N_RealMatrixTranspose" "N_SetRealMatrix"))
 (:c "N_LineFit3dPts" :lisp "line-fit-3d-pts" :file "NL_Geom.cpp" :line 7304 :lines 227 :wave 0 :status :done
   :callees ("N_AllocPt1dArray" "N_AllocReal1dArray" "N_CreatePlanePtNormal" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_IsectPlanePlane" "N_LineFit2dPts" "N_LineImplicitToPtVector" "N_PtFromXYZ" "N_PtToXYZ" "N_RealMatrixSingleValueDecompose" "N_SetRealMatrix" "N_VectorDiff" "N_VectorDot" "N_VectorScale" "N_VectorSum"))
 (:c "N_LineFitPts" :lisp "line-fit-pts" :file "NL_Geom.cpp" :line 7723 :lines 94 :wave 0 :status :done
   :callees ("N_CreateLineStartDirVector" "N_DistPtInfLine" "N_DistSqPtPt2d" "N_FloatOpIsBad" "N_LineFit3dPts" "N_VectorCopy" "N_VectorCreate" "N_VectorDiff" "N_VectorDot" "N_VectorNormalizeRef" "N_VectorReverseInPlace" "N_VectorScale" "N_VectorSum"))
 (:c "N_PtsAreColinear" :lisp "pts-are-colinear" :file "NL_Geom.cpp" :line 9474 :lines 42 :wave 0 :status :done
   :callees ("N_DistPerpPtLineSeg" "N_DistSqPtPt"))
 (:c "N_CalcCircCenterAndRadius" :lisp "calc-circ-center-and-radius" :file "NL_Geom.cpp" :line 9571 :lines 88 :wave 0 :status :done
   :callees ("N_CreateLineStartDirVector" "N_CreatePlanePtNormal" "N_DistPtPt" "N_IsectLineLine" "N_IsectPlanePlane" "N_PtToXYZ" "N_PtsAreColinear" "N_VectorCombine" "N_VectorCopy" "N_VectorCreate" "N_VectorCross" "N_VectorDiff" "N_VectorPerpendicular" "N_VectorScale" "N_VectorSum"))
 (:c "N_CreateCircFrom2dPts" :lisp "create-circ-from-2d-pts" :file "NL_Geom.cpp" :line 9775 :lines 203 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_CalcCircCenterAndRadius" "N_EndNurbs" "N_ErrClear" "N_FloatOpIsBad" "N_GetRealMatrixPtr" "N_InitNurbs" "N_LineFitPts" "N_PtToXYZ" "N_RealMatrixForBackPivot" "N_RealMatrixLstSqSolve" "N_RealMatrixLuDecomposePivot" "N_SetRealMatrix"))
 (:c "N_Combine2CPts" :lisp "combine-2-c-pts" :file "NL_Math.cpp" :line 61 :lines 22 :wave 0 :status :done
   :callees ())
 (:c "N_Combine2Pts" :lisp "combine-2-pts" :file "NL_Math.cpp" :line 182 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_TranslateSum2Pts" :lisp "translate-sum-2-pts" :file "NL_Math.cpp" :line 266 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_TranslateSum2CPts" :lisp "translate-sum-2-c-pts" :file "NL_Math.cpp" :line 311 :lines 27 :wave 0 :status :done
   :callees ())
 (:c "N_Combine4Pts" :lisp "combine-4-pts" :file "NL_Math.cpp" :line 379 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_Combine4CPts" :lisp "combine-4-c-pts" :file "NL_Math.cpp" :line 427 :lines 21 :wave 0 :status :done
   :callees ())
 (:c "N_CopyCPt" :lisp "copy-c-pt" :file "NL_Math.cpp" :line 479 :lines 9 :wave 0 :status :done
   :callees ())
 (:c "N_CopyPt" :lisp "copy-pt" :file "NL_Math.cpp" :line 517 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_CPtToPtEuclid" :lisp "c-pt-to-pt-euclid" :file "NL_Math.cpp" :line 554 :lines 22 :wave 0 :status :done
   :callees ())
 (:c "N_CPtToPtAndW" :lisp "c-pt-to-pt-and-w" :file "NL_Math.cpp" :line 613 :lines 5 :wave 0 :status :done
   :callees ("N_CPtGetW" "N_CPtToPtEuclid"))
 (:c "N_Weight" :lisp "weight" :file "NL_Math.cpp" :line 651 :lines 20 :wave 0 :status :done
   :callees ())
 (:c "N_VectorBlendCPt" :lisp "vector-blend-c-pt" :file "NL_Math.cpp" :line 705 :lines 16 :wave 0 :status :done
   :callees ())
 (:c "N_ScaleCPt" :lisp "scale-c-pt" :file "NL_Math.cpp" :line 753 :lines 16 :wave 0 :status :done
   :callees ())
 (:c "N_ScaleCPtXYZ" :lisp "scale-c-pt-xyz" :file "NL_Math.cpp" :line 802 :lines 10 :wave 0 :status :done
   :callees ())
 (:c "N_VectorBlendPt" :lisp "vector-blend-pt" :file "NL_Math.cpp" :line 847 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_ScalePt" :lisp "scale-pt" :file "NL_Math.cpp" :line 884 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_TransformCPt" :lisp "transform-c-pt" :file "NL_Math.cpp" :line 932 :lines 43 :wave 0 :status :done
   :callees ("N_GetRealMatrixPtr"))
 (:c "N_TranslateCPt" :lisp "translate-c-pt" :file "NL_Math.cpp" :line 1086 :lines 28 :wave 0 :status :done
   :callees ())
 (:c "N_PtToCPt" :lisp "pt-to-c-pt" :file "NL_Math.cpp" :line 1229 :lines 7 :wave 0 :status :done
   :callees ())
 (:c "N_CPtGetW" :lisp "c-pt-get-w" :file "NL_Math.cpp" :line 1267 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_PtToXYZ" :lisp "pt-to-xyz" :file "NL_Math.cpp" :line 1302 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_PtFromXYZ" :lisp "pt-from-xyz" :file "NL_Math.cpp" :line 1339 :lines 6 :wave 0 :status :done :root t
   :callees ())
 (:c "N_CPtToWxWyWz" :lisp "c-pt-to-wx-wy-wz" :file "NL_Math.cpp" :line 1376 :lines 12 :wave 0 :status :done
   :callees ())
 (:c "N_CPtFromWxWyWz" :lisp "c-pt-from-wx-wy-wz" :file "NL_Math.cpp" :line 1420 :lines 12 :wave 0 :status :done
   :callees ())
 (:c "N_CPtSetX" :lisp "c-pt-set-x" :file "NL_Math.cpp" :line 1463 :lines 4 :wave 0 :status :done :root t
   :callees ())
 (:c "N_CPtSetY" :lisp "c-pt-set-y" :file "NL_Math.cpp" :line 1499 :lines 4 :wave 0 :status :done :root t
   :callees ())
 (:c "N_CPtSetZ" :lisp "c-pt-set-z" :file "NL_Math.cpp" :line 1534 :lines 4 :wave 0 :status :done :root t
   :callees ())
 (:c "N_CPtSetW" :lisp "c-pt-set-w" :file "NL_Math.cpp" :line 1570 :lines 4 :wave 0 :status :done :root t
   :callees ())
 (:c "N_Sum2Pts" :lisp "sum-2-pts" :file "NL_Math.cpp" :line 1606 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_Sum2CPts" :lisp "sum-2-c-pts" :file "NL_Math.cpp" :line 1643 :lines 17 :wave 0 :status :done
   :callees ())
 (:c "N_Diff2Pts" :lisp "diff-2-pts" :file "NL_Math.cpp" :line 1741 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_Diff2CPts" :lisp "diff-2-c-pts" :file "NL_Math.cpp" :line 1778 :lines 20 :wave 0 :status :done
   :callees ())
 (:c "N_Dot2Pts" :lisp "dot-2-pts" :file "NL_Math.cpp" :line 1879 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_PtMagnitude" :lisp "pt-magnitude" :file "NL_Math.cpp" :line 2120 :lines 6 :wave 0 :status :done
   :callees ("N_Dot2Pts"))
 (:c "N_CPtRatDeriv" :lisp "c-pt-rat-deriv" :file "NL_Math.cpp" :line 2208 :lines 14 :wave 0 :status :done
   :callees ("N_CPtGetW" "N_CPtToPtEuclid"))
 (:c "N_CPtToPt" :lisp "c-pt-to-pt" :file "NL_Math.cpp" :line 2445 :lines 9 :wave 0 :status :done
   :callees ())
 (:c "N_FloatOpIsBad" :lisp "float-op-is-bad" :file "NL_Math.cpp" :line 2708 :lines 75 :wave 0 :status :done
   :callees ())
 (:c "N_SwapIntegers" :lisp "swap-integers" :file "NL_Math.cpp" :line 2908 :lines 7 :wave 0 :status :vanish
   :callees ())
 (:c "N_SwapReals" :lisp "swap-reals" :file "NL_Math.cpp" :line 2945 :lines 7 :wave 0 :status :done
   :callees ())
 (:c "N_ShellSortReal" :lisp "shell-sort-real" :file "NL_Math.cpp" :line 2983 :lines 20 :wave 0 :status :done
   :callees ())
 (:c "N_SearchInsertReal" :lisp "search-insert-real" :file "NL_Math.cpp" :line 3042 :lines 22 :wave 0 :status :port
   :callees ())
 (:c "N_SwapCPts" :lisp "swap-c-pts" :file "NL_Math.cpp" :line 3394 :lines 7 :wave 0 :status :done
   :callees ("N_CopyCPt"))
 (:c "N_PascalTriRow" :lisp "pascal-tri-row" :file "NL_Math.cpp" :line 3514 :lines 26 :wave 0 :status :done
   :callees ())
 (:c "N_InitRealMatrix" :lisp "init-real-matrix" :file "NL_Math.cpp" :line 3617 :lines 7 :wave 0 :status :vanish
   :callees ())
 (:c "N_CreateRealMatrix" :lisp "create-real-matrix" :file "NL_Math.cpp" :line 3783 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_RealMatrixDefine" :lisp "real-matrix-define" :file "NL_Math.cpp" :line 3966 :lines 7 :wave 0 :status :done
   :callees ())
 (:c "N_GetRealMatrixData" :lisp "get-real-matrix-data" :file "NL_Math.cpp" :line 4139 :lines 8 :wave 0 :status :vanish
   :callees ())
 (:c "N_GetMaxIndexRealMatrix" :lisp "get-max-index-real-matrix" :file "NL_Math.cpp" :line 4298 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_GetRealMatrixPtr" :lisp "get-real-matrix-ptr" :file "NL_Math.cpp" :line 4437 :lines 4 :wave 0 :status :vanish
   :callees ())
 (:c "N_RealMatrixIsNULL" :lisp "real-matrix-is-null" :file "NL_Math.cpp" :line 4707 :lines 15 :wave 0 :status :done
   :callees ("N_GetRealMatrixData"))
 (:c "N_InitRealHash" :lisp "init-real-hash" :file "NL_Math.cpp" :line 4872 :lines 14 :wave 0 :status :port
   :callees ("N_AllocRealPtr1dArray"))
 (:c "N_UpdateRealHash" :lisp "update-real-hash" :file "NL_Math.cpp" :line 4948 :lines 43 :wave 0 :status :port
   :callees ("N_AllocReal1dArray" "N_Realloc1dRealArray"))
 (:c "N_RetrieveRealHash" :lisp "retrieve-real-hash" :file "NL_Math.cpp" :line 5060 :lines 14 :wave 0 :status :port
   :callees ())
 (:c "N_SetRealMatrix" :lisp "set-real-matrix" :file "NL_Math.cpp" :line 5190 :lines 25 :wave 0 :status :done
   :callees ("N_AllocReal2dArray" "N_CreateRealMatrix"))
 (:c "N_CheckMemRealMatrix" :lisp "check-mem-real-matrix" :file "NL_Math.cpp" :line 5493 :lines 18 :wave 0 :status :vanish
   :callees ("N_CheckRealMatStorage" "N_RealMatrixDefine" "N_RealMatrixIsNULL" "N_SetRealMatrix"))
 (:c "N_RealMatrixMultiply" :lisp "real-matrix-multiply" :file "NL_Math.cpp" :line 5714 :lines 145 :wave 0 :status :done
   :callees ("N_CheckMemRealMatrix" "N_GetRealMatrixData" "N_GetRealMatrixPtr"))
 (:c "N_RealMatrixTranspose" :lisp "real-matrix-transpose" :file "NL_Math.cpp" :line 6300 :lines 91 :wave 0 :status :done
   :callees ("N_AllocReal2dArray" "N_CheckMemRealMatrix" "N_CreateRealMatrix" "N_FreeReal2dArray" "N_GetRealMatrixData" "N_GetRealMatrixPtr"))
 (:c "N_RealMatrixMultiplyPtArray" :lisp "real-matrix-multiply-pt-array" :file "NL_Math.cpp" :line 6463 :lines 61 :wave 0 :status :done
   :callees ("N_CopyPt" "N_GetRealMatrixData" "N_VectorBlendPt"))
 (:c "N_RealMatrixMultiplyCPtArray" :lisp "real-matrix-multiply-c-pt-array" :file "NL_Math.cpp" :line 6587 :lines 61 :wave 0 :status :done
   :callees ("N_CopyCPt" "N_GetRealMatrixData" "N_VectorBlendCPt"))
 (:c "N_RealMatrixLuDecompose" :lisp "real-matrix-lu-decompose" :file "NL_Math.cpp" :line 7298 :lines 65 :wave 0 :status :done
   :callees ("N_GetRealMatrixData"))
 (:c "N_RealMatrixLuDecomposePivot" :lisp "real-matrix-lu-decompose-pivot" :file "NL_Math.cpp" :line 7432 :lines 74 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_EndNurbs" "N_GetRealMatrixData" "N_InitNurbs" "N_SwapReals"))
 (:c "N_RealMatrixForBack" :lisp "real-matrix-for-back" :file "NL_Math.cpp" :line 7605 :lines 146 :wave 0 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocPt1dArray" "N_Combine2CPts" "N_Combine2Pts" "N_CopyCPt" "N_CopyPt" "N_Diff2CPts" "N_Diff2Pts" "N_FloatOpIsBad" "N_FreeCPt1dArray" "N_FreePt1dArray" "N_GetRealMatrixData" "N_PtToCPt" "N_VectorBlendCPt" "N_VectorBlendPt"))
 (:c "N_RealMatrixRightForBack" :lisp "real-matrix-right-for-back" :file "NL_Math.cpp" :line 7852 :lines 64 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_EndNurbs" "N_FloatOpIsBad" "N_GetRealMatrixData" "N_InitNurbs"))
 (:c "N_RealMatrixForBackPivot" :lisp "real-matrix-for-back-pivot" :file "NL_Math.cpp" :line 7991 :lines 42 :wave 0 :status :done
   :callees ("N_FloatOpIsBad" "N_GetRealMatrixData"))
 (:c "N_RealMatrixInverse" :lisp "real-matrix-inverse" :file "NL_Math.cpp" :line 8097 :lines 49 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_CheckMemRealMatrix" "N_GetRealMatrixData" "N_GetRealMatrixPtr" "N_RealMatrixLuDecompose" "N_RealMatrixRightForBack"))
 (:c "N_RealMatrixLstSqSolve" :lisp "real-matrix-lst-sq-solve" :file "NL_Math.cpp" :line 8337 :lines 43 :wave 0 :status :done
   :callees ("N_AllocInt1dArray" "N_EndNurbs" "N_GetRealMatrixData" "N_InitNurbs" "N_InitRealMatrix" "N_RealMatrixForBackPivot" "N_RealMatrixLuDecomposePivot" "N_RealMatrixMultiply" "N_RealMatrixTranspose"))
 (:c "N_RealMatrixSingleValueDecompose" :lisp "real-matrix-single-value-decompose" :file "NL_Math.cpp" :line 8707 :lines 254 :wave 0 :status :done
   :callees ("N_AllocReal1dArray" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_FuncFindMinima" :lisp "func-find-minima" :file "NL_Math.cpp" :line 9416 :lines 110 :wave 0 :status :done
   :callees ())
 (:c "N_CrvApproxSpiral" :lisp "crv-approx-spiral" :file "NL_Spiral.cpp" :line 88 :lines 121 :wave 5 :status :done :root t
   :callees ("N_AllocPt1dArray" "N_CrvEval" "N_CrvGetKnots" "N_EndNurbs" "N_FitCubicSplineInterp" "N_InitNurbs" "ST_GetPointAtAngle" "ST_GetPointAtAngleD"))
 (:c "ST_GetSpiralAngle" :lisp "spiral--st-get-spiral-angle" :file "NL_Spiral.cpp" :line 459 :lines 4 :wave 5 :status :done :static t
   :callees ())
 (:c "ST_GetPointAtAngle" :lisp "spiral--st-get-point-at-angle" :file "NL_Spiral.cpp" :line 468 :lines 12 :wave 5 :status :done :static t
   :callees ("ST_GetSpiralAngle"))
 (:c "ST_GetPointAtAngleD" :lisp "spiral--st-get-point-at-angle-d" :file "NL_Spiral.cpp" :line 485 :lines 13 :wave 5 :status :done :static t
   :callees ("ST_GetSpiralAngle"))
 (:c "N_SrfEvalPtGrid" :lisp "srf-eval-pt-grid" :file "NL_SrfAdv.cpp" :line 94 :lines 96 :wave 6 :status :done :root t
   :callees ("N_AllocCPt1dArray" "N_AllocInt1dArray" "N_AllocReal2dArray" "N_BasisEvalArray" "N_CPtToPtEuclid" "N_CopyCPt" "N_EndNurbs" "N_InitNurbs" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_VectorBlendCPt"))
 (:c "N_SrfExtendToCrv" :lisp "srf-extend-to-crv" :file "NL_SrfAdv.cpp" :line 2711 :lines 643 :wave 6 :status :done :root t :fenced t
   :callees ("N_AllocCPt1dArray" "N_AllocCPt2dArray" "N_AllocSrfArrays" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_CrvCopy" "N_CrvEval" "N_CrvExtendToPt" "N_CrvGetCPts" "N_CrvGetParamBounds" "N_CrvInitArrays" "N_CrvMakeCompatibleWithSrf" "N_CrvSizeArrays" "N_EndNurbs" "N_FreeCrv" "N_FreeSrf" "N_InitNurbs" "N_IsSrfRat" "N_ScaleCPt" "N_SrfCopy" "N_SrfExtractIsoCrv" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetParameterBounds" "N_SrfInitArrays" "N_SrfInsertKnot" "N_SrfRemoveKnotConditional" "N_SrfSizeArrays" "N_VectorBlendCPt"))
 (:c "N_SrfModifyBoundaryCrv" :lisp "srf-modify-boundary-crv" :file "NL_SrfAdv.cpp" :line 3782 :lines 191 :wave 6 :status :done :root t
   :callees ("N_AllocSrfArrays" "N_CPtFromWxWyWz" "N_CopyCPt" "N_CreateRectangle" "N_CrvCopy" "N_CrvGetCPtsDegreeAndKnots" "N_CrvInitArrays" "N_CrvMakeCompatibleWithSrf" "N_Diff2CPts" "N_EndNurbs" "N_FreeSrf" "N_InitNurbs" "N_IsCrvRat" "N_MakeSrfsCompatibleUV" "N_SrfCopy" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetParameterBounds" "N_SrfInitArrays" "N_SrfReparamToInterval" "N_SrfSizeArrays" "N_Sum2CPts"))
 (:c "N_SrfEvalPtNormalDeriv" :lisp "srf-eval-pt-normal-deriv" :file "NL_SrfAdv.cpp" :line 4455 :lines 102 :wave 6 :status :done
   :callees ("N_AllocPt2dArray" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_SrfDerivs" "N_VectorCross" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude" "N_VectorScale" "N_VectorSum"))
 (:c "N_SrfOffsetGetMaxDeriv" :lisp "srf-offset-get-max-deriv" :file "NL_SrfAdv.cpp" :line 4642 :lines 90 :wave 6 :status :done
   :callees ("N_AllocPt2dArray" "N_AllocReal1dArray" "N_EndNurbs" "N_ErrClear" "N_InitNurbs" "N_PtMagnitude" "N_SrfDecomposeToBez" "N_SrfDerivs" "N_SrfEvalPtNormalDeriv" "N_SrfGetDegrees" "N_SrfGetKnots"))
 (:c "N_SrfGetGridPtsNormals" :lisp "srf-get-grid-pts-normals" :file "NL_SrfAdv.cpp" :line 4825 :lines 153 :wave 6 :status :done
   :callees ("N_AllocInt1dArray" "N_AllocPt2dArray" "N_AllocReal2dArray" "N_AllocReal3dArray" "N_BasisDerivsArray" "N_CPtGetW" "N_CPtToPt" "N_CopyPt" "N_Diff2Pts" "N_EndNurbs" "N_InitNurbs" "N_IsSrfRat" "N_ScalePt" "N_SrfEvalPtNormalAtPole" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_VectorBlendPt" "N_VectorCross" "N_VectorNormalize"))
 (:c "N_ReparmCrvsIsectPt" :lisp "reparm-crvs-isect-pt" :file "NL_SrfAdv.cpp" :line 7521 :lines 108 :wave 6 :status :done :root t
   :callees ("N_AllocCFuncArrays" "N_AllocReal1dArray" "N_CFuncDerivs" "N_CFuncGetKnots" "N_CFuncSetSizeIndices" "N_CrvsMakeCompatible" "N_EndNurbs" "N_FitFuncInterpGivenParams" "N_InitNurbs" "N_SrfReparamFunc"))
 (:c "N_CreateSkinSrf" :lisp "create-skin-srf" :file "NL_SrfAdv.cpp" :line 7899 :lines 254 :wave 6 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocCPt2dArray" "N_AllocCrvArrays" "N_AllocKnotVectorAndArray" "N_AllocReal1dArray" "N_CopyCPt" "N_CrvDerivsAtKnot" "N_CrvGetCPts" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvSetSizeIndices" "N_CrvsAreCoincident" "N_CrvsMakeCompatibleKnotTol" "N_EndNurbs" "N_FitCalcKnotVectorEndDerivs" "N_FitCalcSrfParamValues" "N_FitCrvCalcKnotVector" "N_FitCrvInterpGivenParams" "N_FitCrvKnotsAndDerivs" "N_InitNurbs" "N_SrfGetCPtsAndKnots" "N_SrfGetKnotVectors" "N_SrfSizeArrays" "N_SrfWeightsAreValid"))
 (:c "N_CreateSkinSrfParams" :lisp "create-skin-srf-params" :file "NL_SrfAdv.cpp" :line 8800 :lines 227 :wave 6 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocCPt2dArray" "N_AllocCrvArrays" "N_AllocKnotVectorAndArray" "N_AllocReal1dArray" "N_CopyCPt" "N_CrvDerivsAtKnot" "N_CrvGetCPts" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvSetSizeIndices" "N_CrvsAreCoincident" "N_CrvsMakeCompatible" "N_EndNurbs" "N_FitCalcSrfParamValues" "N_FitCrvCalcKnotVector" "N_FitCrvInterpGivenParams" "N_FitCrvKnotsAndDerivs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays" "N_SrfWeightsAreValid"))
 (:c "N_CreateSkinSpine" :lisp "create-skin-spine" :file "NL_SrfAdv.cpp" :line 10401 :lines 479 :wave 6 :status :done :root t
   :callees ("N_AllocCFuncArrays" "N_AllocCPt1dArray" "N_AllocCPt2dArray" "N_AllocCrvArrays" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_CFuncDerivs" "N_CPtGetW" "N_CPtRatDeriv" "N_CopyCPt" "N_CopyPt" "N_CreateBVectors" "N_CreateSkinSrf" "N_CreateTransformMatrixFromAxes" "N_CrvDerivs" "N_CrvEval" "N_CrvEvalTangent" "N_CrvGetCPts" "N_CrvGetCPtsAndKnots" "N_CrvGetCPtsDegreeAndKnots" "N_CrvTransform" "N_CrvsAreCoincident" "N_CrvsMakeCompatible" "N_DistCptCpt" "N_DistPtPt" "N_EndNurbs" "N_FitCalcKnotVectorEndDerivs" "N_FitCalcKnotsDerivs" "N_FitCalcSrfParamValues" "N_FitCrvFirstDerivAndKnots" "N_FitCrvKnotsAndDerivs" "N_FitFuncInterpGivenParams" "N_InitNurbs" "N_IsCrvRat" "N_PtToCPt" "N_ScalePt" "N_SetRealMatrix" "N_SrfGetCPtsAndKnots" "N_SrfGetKnotVectors" "N_SrfSizeArrays" "N_SrfWeightsAreValid" "N_VectorCross"))
 (:c "N_CreateBVectors" :lisp "create-b-vectors" :file "NL_SrfAdv.cpp" :line 11270 :lines 50 :wave 6 :status :done
   :callees ("N_AllocPt1dArray" "N_Combine2Pts" "N_CopyPt" "N_CrvEvalTangent" "N_CrvIsClosed" "N_DistPtPt" "N_EndNurbs" "N_InitNurbs" "N_VectorCopy" "N_VectorDot" "N_VectorNormalizeRef"))
 (:c "N_CreateGordonSrf" :lisp "create-gordon-srf" :file "NL_SrfAdv.cpp" :line 11775 :lines 120 :wave 6 :status :done :root t
   :callees ("N_AllocArraySrfPtrs" "N_AllocKnotVectorAndArray" "N_AllocPt2dArray" "N_Combine2Pts" "N_CreateSkinSrf" "N_CreateSkinSrfParams" "N_CrvEval" "N_CrvsAreCombatible" "N_FitCrvCalcKnotVector" "N_FitSrfToPts" "N_FitSrfToPtsKnots" "N_InitNurbs" "N_IsCrvRat" "N_MakeSrfsCompatibleUV" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfInitArrays" "N_SrfSizeArrays" "N_TranslateSum2CPts"))
 (:c "N_CreateCoonsBoundaryCrvs" :lisp "create-coons-boundary-crvs" :file "NL_SrfAdv.cpp" :line 17715 :lines 321 :wave 6 :status :done :root t
   :callees ("N_AllocArrayCrvPtrs" "N_AllocArraySrfPtrs" "N_AllocSrfArrays" "N_Combine2CPts" "N_Combine4CPts" "N_CopyCPt" "N_CrvDerivs" "N_CrvEval" "N_CrvGetCPts" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnots" "N_CrvGetParamBounds" "N_CrvsAreCombatible" "N_CrvsMakeCompatible" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat" "N_MakeSrfsCompatibleUV" "N_PtToCPt" "N_SrfGetCPts" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfSizeArrays" "N_TranslateSum2CPts" "N_VectorCombine" "N_VectorDiff" "N_VectorMagnitude"))
 (:c "N_ApproxSubSrfWithSrf" :lisp "approx-sub-srf-with-srf" :file "NL_SrfApprox.cpp" :line 1243 :lines 519 :wave 6 :status :done :root t
   :callees ("N_AllocArrayCrvPtrsAndData" "N_AllocCrvFuncArray" "N_AllocPt2dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisFindIndexNodeArray" "N_CFuncEval" "N_CopyPt" "N_CrvClosestPt" "N_CrvCopy" "N_CrvElevateDegree" "N_CrvEval" "N_CrvFuncCopy" "N_CrvGetArraySizes" "N_CrvGetDegree" "N_CrvGetKnotVector" "N_CrvGetParamBounds" "N_CrvRatToNonRat" "N_CrvReparam" "N_CrvsMakeCompatible" "N_DistPtPt" "N_EndNurbs" "N_FitFuncInterpGivenParams" "N_FitSrfFuncInterp" "N_FitSrfToPtsAndBoundary" "N_FuncFindMinima" "N_GetClosestPtOnSrf" "N_InitNurbs" "N_PtFromXYZ" "N_PtToXYZ" "N_SrfEvalPt" "N_SrfGetParameterBounds" "N_SrfInitArrays" "ST_GetDist" "ST_GetPoint"))
 (:c "ST_GetDist" :lisp "srf-approx--st-get-dist" :file "NL_SrfApprox.cpp" :line 2612 :lines 23 :wave 6 :status :done :static t
   :callees ("N_CrvEval" "N_DistPtPt" "N_PtToXYZ" "N_SrfEvalPt"))
 (:c "ST_GetPoint" :lisp "srf-approx--st-get-point" :file "NL_SrfApprox.cpp" :line 2652 :lines 54 :wave 6 :status :done :static t
   :callees ("N_CFuncEval" "N_CrvEval" "N_PtFromXYZ" "N_PtToXYZ"))
 (:c "N_SrfWeightsAreValid" :lisp "srf-weights-are-valid" :file "NL_SrfBasic.cpp" :line 132 :lines 27 :wave 0 :status :done
   :callees ("N_CPtGetW" "N_ErrSet" "N_IsSrfRat" "N_SrfGetCPts"))
 (:c "N_SrfIsSized" :lisp "srf-is-sized" :file "NL_SrfBasic.cpp" :line 285 :lines 12 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_SrfGetArraySizes"))
 (:c "N_AllocSrf" :lisp "alloc-srf" :file "NL_SrfBasic.cpp" :line 414 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocSrfAndArrays" :lisp "alloc-srf-and-arrays" :file "NL_SrfBasic.cpp" :line 489 :lines 20 :wave 0 :status :vanish
   :callees ("N_AllocCNetAndArrays" "N_AllocKnotVectorAndArray" "N_AllocSrf" "N_SrfFromCNetAndKnotVectors"))
 (:c "N_AllocArraySrfPtrs" :lisp "alloc-array-srf-ptrs" :file "NL_SrfBasic.cpp" :line 559 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_Alloc2dArraySrfPtrsParameters" :lisp "alloc-2d-array-srf-ptrs-parameters" :file "NL_SrfBasic.cpp" :line 751 :lines 72 :wave 0 :status :vanish
   :callees ("N_AllocSrfAndArrays" "N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocArraySrfPtrsInit" :lisp "alloc-array-srf-ptrs-init" :file "NL_SrfBasic.cpp" :line 888 :lines 33 :wave 0 :status :vanish :root t
   :callees ("N_AllocSrf" "N_ErrSet" "N_Free" "N_Malloc" "N_SrfInitArrays"))
 (:c "N_FreeSrf" :lisp "free-srf" :file "NL_SrfBasic.cpp" :line 971 :lines 20 :wave 0 :status :vanish
   :callees ("N_FreeCNet" "N_FreeCPt2dArray" "N_FreeKnotVector" "N_FreeReal1dArray" "N_SrfGetCPtsAndKnots" "N_SrfGetNetAndKnotVectors"))
 (:c "N_SrfFromCNetAndKnotVectors" :lisp "srf-from-c-net-and-knot-vectors" :file "NL_SrfBasic.cpp" :line 1108 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_SrfFromCPtsAndKnots" :lisp "srf-from-c-pts-and-knots" :file "NL_SrfBasic.cpp" :line 1160 :lines 19 :wave 0 :status :done :root t
   :callees ("N_AllocCNet" "N_AllocKnotVector" "N_CNetFromCPts" "N_KnotVectorFromRealArray" "N_SrfFromCNetAndKnotVectors"))
 (:c "N_SrfFromCPtCoordsAndKnots" :lisp "srf-from-c-pt-coords-and-knots" :file "NL_SrfBasic.cpp" :line 1237 :lines 23 :wave 0 :status :done :root t
   :callees ("N_AllocCNet" "N_AllocKnotVector" "N_CNetFromCPts" "N_KnotVectorFromRealArray" "N_SrfFromCNetAndKnotVectors" "N_XYZTo2dCPtArray"))
 (:c "N_SrfGetNetAndKnotVectors" :lisp "srf-get-net-and-knot-vectors" :file "NL_SrfBasic.cpp" :line 1409 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_SrfGetCPtsDegreesAndKnots" :lisp "srf-get-c-pts-degrees-and-knots" :file "NL_SrfBasic.cpp" :line 1453 :lines 22 :wave 0 :status :vanish
   :callees ())
 (:c "N_SrfGetCPts" :lisp "srf-get-c-pts" :file "NL_SrfBasic.cpp" :line 1506 :lines 6 :wave 0 :status :vanish
   :callees ())
 (:c "N_SrfGetArraySizes" :lisp "srf-get-array-sizes" :file "NL_SrfBasic.cpp" :line 1541 :lines 7 :wave 0 :status :vanish
   :callees ())
 (:c "N_SrfGetCPtsAndKnots" :lisp "srf-get-c-pts-and-knots" :file "NL_SrfBasic.cpp" :line 1579 :lines 10 :wave 0 :status :vanish
   :callees ())
 (:c "N_SrfGetCPtsKnotVectorAndKnots" :lisp "srf-get-c-pts-knot-vector-and-knots" :file "NL_SrfBasic.cpp" :line 1620 :lines 8 :wave 0 :status :done
   :callees ())
 (:c "N_SrfSetCPtsAndKnots" :lisp "srf-set-c-pts-and-knots" :file "NL_SrfBasic.cpp" :line 1659 :lines 6 :wave 0 :status :done
   :callees ())
 (:c "N_SrfGetKnots" :lisp "srf-get-knots" :file "NL_SrfBasic.cpp" :line 1695 :lines 7 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_SrfGetParameterBounds" :lisp "srf-get-parameter-bounds" :file "NL_SrfBasic.cpp" :line 1733 :lines 7 :wave 0 :status :done
   :callees ())
 (:c "N_SrfGetKnotVectors" :lisp "srf-get-knot-vectors" :file "NL_SrfBasic.cpp" :line 1768 :lines 5 :wave 0 :status :vanish
   :callees ())
 (:c "N_SrfGetDegrees" :lisp "srf-get-degrees" :file "NL_SrfBasic.cpp" :line 1801 :lines 5 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_SrfCompress" :lisp "srf-compress" :file "NL_SrfBasic.cpp" :line 1842 :lines 35 :wave 0 :status :done
   :callees ("N_AllocCPt2dArray" "N_AllocReal1dArray" "N_CopyCPt" "N_FreeCPt2dArray" "N_FreeReal1dArray" "N_SrfGetArraySizes" "N_SrfGetCPtsAndKnots" "N_SrfSetCPtsAndKnots"))
 (:c "N_ClampSrfAtParams" :lisp "clamp-srf-at-params" :file "NL_SrfBasic.cpp" :line 1932 :lines 14 :wave 0 :status :port
   :callees ("N_SrfGetKnots"))
 (:c "N_SrfSetSizeIndices" :lisp "srf-set-size-indices" :file "NL_SrfBasic.cpp" :line 1983 :lines 9 :wave 0 :status :vanish
   :callees ())
 (:c "N_SwapUV" :lisp "swap-uv" :file "NL_SrfBasic.cpp" :line 2025 :lines 24 :wave 0 :status :done :root t
   :callees ("N_AllocCPt2dArray" "N_CNetFromCPts" "N_CopyCPt" "N_FreeCPt2dArray" "N_SrfFromCNetAndKnotVectors" "N_SrfGetCPts" "N_SrfGetNetAndKnotVectors"))
 (:c "N_SrfInitArrays" :lisp "srf-init-arrays" :file "NL_SrfBasic.cpp" :line 2100 :lines 8 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_SrfAreArraysNULL" :lisp "srf-are-arrays-null" :file "NL_SrfBasic.cpp" :line 2137 :lines 15 :wave 0 :status :vanish
   :callees ("N_SrfGetNetAndKnotVectors"))
 (:c "N_IsSrfRat" :lisp "is-srf-rat" :file "NL_SrfBasic.cpp" :line 2188 :lines 16 :wave 0 :status :done :root t
   :callees ("N_CPtGetW" "N_SrfGetCPts"))
 (:c "N_SrfIsDegen" :lisp "srf-is-degen" :file "NL_SrfBasic.cpp" :line 2240 :lines 39 :wave 0 :status :port
   :callees ("N_CPtToPtEuclid" "N_CopyPt" "N_DistPtPt" "N_ScalePt" "N_SrfGetCPts" "N_Sum2Pts"))
 (:c "N_SrfIsClosed" :lisp "srf-is-closed" :file "NL_SrfBasic.cpp" :line 2323 :lines 41 :wave 0 :status :done
   :callees ("N_CPtToPtEuclid" "N_DistPtPt" "N_SrfGetCPts"))
 (:c "N_1dPtSetIsClosed" :lisp "1d-pt-set-is-closed" :file "NL_SrfBasic.cpp" :line 2420 :lines 25 :wave 0 :status :done
   :callees ("N_DistCptCptHomo" "N_DistPtPt"))
 (:c "N_SrfIsSingular" :lisp "srf-is-singular" :file "NL_SrfBasic.cpp" :line 2858 :lines 126 :wave 0 :status :port
   :callees ("N_CPtToPtEuclid" "N_DistPtPt" "N_ScalePt" "N_SrfGetCPts" "N_Sum2Pts"))
 (:c "N_SrfIsFlat" :lisp "srf-is-flat" :file "NL_SrfBasic.cpp" :line 3069 :lines 292 :wave 0 :status :port
   :callees ("N_AllocPt2dArray" "N_AllocReal2dArray" "N_CPtToPtEuclid" "N_Combine4Pts" "N_CrvInitArrays" "N_CrvIsLine" "N_DistCptCpt" "N_DistPtLineSeg" "N_DistPtPt" "N_DistSignedPtPlane" "N_EndNurbs" "N_InitNurbs" "N_PtToXYZ" "N_SrfEvalPtPtDerivNormalFast" "N_SrfExtractIsoCrv" "N_SrfGetCPts" "N_SrfGetParameterBounds" "N_SrfIsClosed" "N_SrfIsDegen" "N_SrfIsSingular" "N_VectorCross" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude" "N_VectorNormalizeRef" "N_VectorScale"))
 (:c "N_AllocSrfArrays" :lisp "alloc-srf-arrays" :file "NL_SrfBasic.cpp" :line 3518 :lines 24 :wave 0 :status :vanish
   :callees ("N_AllocCNetAndArrays" "N_AllocKnotVectorAndArray" "N_SrfFromCNetAndKnotVectors"))
 (:c "N_SrfSizeArrays" :lisp "srf-size-arrays" :file "NL_SrfBasic.cpp" :line 3601 :lines 31 :wave 0 :status :vanish
   :callees ("N_AllocSrfArrays" "N_SrfAreArraysNULL" "N_SrfIsSized" "N_SrfSetSizeIndices"))
 (:c "N_SrfCopy" :lisp "srf-copy" :file "NL_SrfBasic.cpp" :line 3684 :lines 24 :wave 0 :status :done
   :callees ("N_CopyCPt" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfSizeArrays"))
 (:c "N_SrfGetENet" :lisp "srf-get-e-net" :file "NL_SrfBasic.cpp" :line 4697 :lines 34 :wave 0 :status :done
   :callees ("N_AllocPt2dArray" "N_CPtToPtEuclid" "N_ENetFromPts" "N_ErrSet" "N_SrfGetCPts"))
 (:c "N_SrfNonRatToRat" :lisp "srf-non-rat-to-rat" :file "NL_SrfBasic.cpp" :line 4971 :lines 13 :wave 0 :status :done
   :callees ("N_CPtSetW" "N_SrfGetCPts"))
 (:c "N_MakeSrfsRatCompatible" :lisp "make-srfs-rat-compatible" :file "NL_SrfBasic.cpp" :line 5026 :lines 24 :wave 0 :status :done
   :callees ("N_IsSrfRat" "N_SrfNonRatToRat"))
 (:c "N_SrfReparamToInterval" :lisp "srf-reparam-to-interval" :file "NL_SrfBasic.cpp" :line 5095 :lines 37 :wave 0 :status :done
   :callees ("N_RectangleGetData" "N_SrfGetDegrees" "N_SrfGetKnots"))
 (:c "N_SrfIsFlatCheap" :lisp "srf-is-flat-cheap" :file "NL_SrfBasic.cpp" :line 5646 :lines 69 :wave 0 :status :done
   :callees ("N_CPtToPtEuclid" "N_Combine4Pts" "N_DistSignedPtPlane" "N_PlanePtNormalToImplicit" "N_SrfGetCPts" "N_SrfIsClosed" "N_VectorCombine" "N_VectorCross" "N_VectorDiff" "N_VectorMagnitude"))
 (:c "N_CreateSrfCornerPts" :lisp "create-srf-corner-pts" :file "NL_SrfCommon.cpp" :line 55 :lines 25 :wave 4 :status :done :root t
   :callees ("N_PtToCPt" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays"))
 (:c "N_CreateSrfExtrudeCrv" :lisp "create-srf-extrude-crv" :file "NL_SrfCommon.cpp" :line 148 :lines 83 :wave 4 :status :done
   :callees ("N_CPtToPtAndW" "N_CopyCPt" "N_Crv2dTo3d" "N_CrvGetCPtsDegreeAndKnots" "N_CrvIs3d" "N_IsCrvRat" "N_PtToCPt" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays" "N_VectorNormalizeRef" "N_VectorPtAlongVector" "N_Weight"))
 (:c "N_CreateRuledSrfFromBoundaryCrvs" :lisp "create-ruled-srf-from-boundary-crvs" :file "NL_SrfCommon.cpp" :line 1046 :lines 186 :wave 4 :status :done :root t :fenced t
   :callees ("N_AllocKnotVectorAndArray" "N_BasisGetSpanCount" "N_CFuncInitArrays" "N_CopyCPt" "N_CreateCrvFromNumAndDenom" "N_CrvCopy" "N_CrvElevateDegree" "N_CrvFuncMultiplyCrv" "N_CrvFuncMultiplyCrvFunc" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetKnotVector" "N_CrvGetNumAndDenom" "N_CrvInitArrays" "N_CrvNonRatToRat" "N_CrvRatToNonRat" "N_CrvRefine" "N_EndNurbs" "N_GetCompatibleKnotArrayMult" "N_InitNurbs" "N_IsCrvRat" "N_KnotVectorGetKnots" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays"))
 (:c "N_CreateRevolvedSrf" :lisp "create-revolved-srf" :file "NL_SrfCommon.cpp" :line 1398 :lines 74 :wave 4 :status :done :root t
   :callees ("N_AllocCrvArrays" "N_CPtGetW" "N_CPtToPtEuclid" "N_CalcCircWeights" "N_CreateCircArc" "N_CreateLineStartDirVector" "N_Crv2dTo3d" "N_CrvGetCPts" "N_CrvGetCPtsDegreeAndKnots" "N_CrvGetConicData" "N_CrvGetKnots" "N_CrvIs3d" "N_CrvNonRatToRat" "N_EndNurbs" "N_InitNurbs" "N_IsCrvRat" "N_ProjectPtLine" "N_ScaleCPt" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays" "N_VectorCross" "N_VectorDiff" "N_VectorMagnitude" "N_Weight"))
 (:c "N_CreateSphere" :lisp "create-sphere" :file "NL_SrfCommon.cpp" :line 1594 :lines 22 :wave 4 :status :done :root t :fenced t
   :callees ("N_CreateCircArc" "N_CreateRevolvedSrf" "N_CrvInitArrays" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_FitCalcSrfParamValues" :lisp "fit-calc-srf-param-values" :file "NL_SrfFit.cpp" :line 101 :lines 124 :wave 5 :status :done
   :callees ("N_AllocReal1dArray" "N_DistCptCptHomo" "N_DistPtPt" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_FitSrfToPts" :lisp "fit-srf-to-pts" :file "NL_SrfFit.cpp" :line 341 :lines 127 :wave 5 :status :done :root t
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_BasisEval" "N_CopyCPt" "N_EndNurbs" "N_FitCalcSrfParamValues" "N_FitCrvCPtsFromSrfData" "N_FitCrvCalcKnotVector" "N_GetRealMatrixPtr" "N_InitNurbs" "N_PtToCPt" "N_RealMatrixLuDecompose" "N_SetRealMatrix" "N_SrfGetCPtsKnotVectorAndKnots" "N_SrfSizeArrays"))
 (:c "N_FitSrfToPtsKnots" :lisp "fit-srf-to-pts-knots" :file "NL_SrfFit.cpp" :line 599 :lines 149 :wave 5 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_BasisEval" "N_CopyCPt" "N_EndNurbs" "N_FitCrvCPtsFromSrfData" "N_GetRealMatrixPtr" "N_InitNurbs" "N_KnotVectorGetKnots" "N_PtToCPt" "N_RealMatrixLuDecompose" "N_SetRealMatrix" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays"))
 (:c "N_FitSrfLstSqApprox" :lisp "fit-srf-lst-sq-approx" :file "NL_SrfFit.cpp" :line 883 :lines 384 :wave 5 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocCPt2dArray" "N_AllocInt1dArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisEval" "N_BasisIEval" "N_CopyCPt" "N_CopyPt" "N_CreateRealMatrix" "N_EndNurbs" "N_FitCalcKnotVectorCrvApprox" "N_FitCalcSrfParamValues" "N_InitNurbs" "N_PtToCPt" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_SrfGetCPtsAndKnots" "N_SrfGetKnotVectors" "N_SrfSizeArrays" "N_TranslateSum2CPts" "N_TranslateSum2Pts" "N_VectorBlendCPt" "N_VectorBlendPt"))
 (:c "N_FitPtsNormals" :lisp "fit-pts-normals" :file "NL_SrfFit.cpp" :line 3513 :lines 508 :wave 5 :status :done :root t
   :callees ("N_AllocPt2dArray" "N_AllocReal1dArray" "N_Combine4Pts" "N_CopyPt" "N_CreatePlanePtNormal" "N_ENetFromPts" "N_ENetIsClosed" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_ProjectPtPlane" "N_PtToCPt" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays" "N_VectorCombine" "N_VectorCopy" "N_VectorCross" "N_VectorDiff" "N_VectorMagnitude" "N_VectorNormalizeRef" "N_VectorPtAlongVector" "N_VectorScale"))
 (:c "N_FitSrfToPtsAndBoundary" :lisp "fit-srf-to-pts-and-boundary" :file "NL_SrfFit.cpp" :line 5576 :lines 650 :wave 5 :status :done
   :callees ("N_AllocArrayCrvPtrs" "N_AllocKnotVectorAndArray" "N_AllocReal1dArray" "N_CPtToPtEuclid" "N_CrvClosestPt" "N_CrvGetCPts" "N_CrvGetCPtsDegreeAndKnots" "N_CrvInitArrays" "N_CrvRefine" "N_EndNurbs" "N_FitCalcKnotVectorCrvApprox" "N_FitSrfLstSqKnots" "N_InitNurbs" "N_KnotVectorFromRealArray" "N_KnotVectorGetKnots" "N_PtToCPt" "N_SetKnotIndex" "N_SrfGetCPts" "N_SrfGetKnots" "N_SrfInsertKnots"))
 (:c "N_FitSrfApproxTol" :lisp "fit-srf-approx-tol" :file "NL_SrfFit.cpp" :line 6545 :lines 228 :wave 5 :status :done :root t
   :callees ("N_AllocPt2dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisGetSpanCount" "N_DistPtPt" "N_EndNurbs" "N_FitCalcKnotVectorCrvApprox" "N_FitCalcSrfParamValues" "N_FitCrvCalcKnotVector" "N_FitSrfApproxRemoveKnots" "N_FitSrfLstSqKnots" "N_FitSrfToPtsKnots" "N_GetClosestPtOnSrf" "N_InitNurbs" "N_Pts2dGetMaxExtent" "N_SrfAreArraysNULL" "N_SrfCompress" "N_SrfGetArraySizes" "N_SrfGetKnotVectors" "N_SrfGetKnots" "N_SrfSetSizeIndices" "N_SrfSizeArrays"))
 (:c "N_FitSrfLstSqKnots" :lisp "fit-srf-lst-sq-knots" :file "NL_SrfFit.cpp" :line 8734 :lines 478 :wave 5 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocCPt2dArray" "N_AllocInt1dArray" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisEval" "N_BasisIEval" "N_CopyCPt" "N_CopyPt" "N_CreateRealMatrix" "N_EndNurbs" "N_InitNurbs" "N_KnotVectorGetKnots" "N_PtToCPt" "N_RealMatrixForBack" "N_RealMatrixLuDecompose" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays" "N_TranslateSum2CPts" "N_TranslateSum2Pts" "N_VectorBlendCPt" "N_VectorBlendPt"))
 (:c "N_FitSrfFuncInterp" :lisp "fit-srf-func-interp" :file "NL_SrfFit.cpp" :line 9438 :lines 125 :wave 5 :status :done
   :callees ("N_AllocReal1dArray" "N_BasisEval" "N_EndNurbs" "N_FitCrvCalcKnotVector" "N_GetRealMatrixPtr" "N_InitNurbs" "N_RealMatrixLuDecompose" "N_RealMatrixRightForBack" "N_SFuncGetKnotVectors" "N_SFuncSizeArrays" "N_SetRealMatrix" "N_SrfFuncCntrlVal"))
 (:c "N_FitSrfInterpBoundary" :lisp "fit-srf-interp-boundary" :file "NL_SrfFit.cpp" :line 9690 :lines 218 :wave 5 :status :done :root t
   :callees ("N_AllocArrayCrvPtrs" "N_AllocCrv" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_CPtToPtEuclid" "N_CopyCPt" "N_CrvCopy" "N_CrvGetArraySizes" "N_CrvGetCPts" "N_CrvGetDegree" "N_CrvGetKnotVector" "N_CrvGetKnots" "N_CrvInitArrays" "N_CrvProjectPts" "N_CrvsAreCombatible" "N_CrvsMakeCompatible" "N_EndNurbs" "N_FitCrvKnotsAndTangents" "N_InitNurbs" "N_IsCrvRat" "N_KnotVectorGetKnots" "N_KnotsCopy" "N_SrfGetCPtsAndKnots" "N_SrfSizeArrays" "N_VectorCopy"))
 (:c "N_FitSrfInterpBicubic" :lisp "fit-srf-interp-bicubic" :file "NL_SrfFit.cpp" :line 12645 :lines 10 :wave 5 :status :done :root t
   :callees ("N_FitPtsNormals"))
 (:c "N_FitSrfApproxRemoveKnots" :lisp "fit-srf-approx-remove-knots" :file "NL_SrfFit.cpp" :line 14164 :lines 626 :wave 5 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocInt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisIEval" "N_CFuncEval" "N_CFuncFromKnots" "N_CFuncSetSizeIndices" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_FitSrfRemovalBoundary" "N_InitNurbs" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_SrfSetSizeIndices"))
 (:c "N_FitSrfRemovalBoundary" :lisp "fit-srf-removal-boundary" :file "NL_SrfFit.cpp" :line 15087 :lines 133 :wave 5 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_Combine2CPts" "N_CopyCPt" "N_DistCptCptHomo" "N_EndNurbs" "N_InitNurbs" "N_SrfGetCPtsDegreesAndKnots"))
 (:c "N_GetClosestPtOnSrf" :lisp "get-closest-pt-on-srf" :file "NL_SrfGeom.cpp" :line 146 :lines 269 :wave 2 :status :done
   :callees ("N_Combine2Pts" "N_CopyPt" "N_DistSqPtPt" "N_FloatOpIsBad" "N_PtMagnitude" "N_SrfDerivs" "N_SrfEvalPt" "N_SrfGetKnots" "N_SrfIsClosed" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude" "N_VectorScale"))
 (:c "N_SrfGetAverageLen" :lisp "srf-get-average-len" :file "NL_SrfGeom.cpp" :line 1316 :lines 50 :wave 2 :status :done
   :callees ("N_AllocPt1dArray" "N_CPtToPtEuclid" "N_DistPtPt" "N_EndNurbs" "N_InitNurbs" "N_SrfGetCPts"))
 (:c "N_ApproxSrfWithQuadSrf" :lisp "approx-srf-with-quad-srf" :file "NL_SrfGeom.cpp" :line 1453 :lines 203 :wave 2 :status :done :root t
   :callees ("N_AllocPt2dArray" "N_AllocReal1dArray" "N_AllocSrfArrays" "N_EndNurbs" "N_InitNurbs" "N_Realloc1dRealArray" "N_SrfEvalPt" "N_SrfExtractPatch" "N_SrfGetArraySizes" "N_SrfGetDegrees" "N_SrfGetKnots" "N_SrfIsFlatCheap" "N_SrfMaxDiagDistBBox" "N_SrfSetSizeIndices"))
 (:c "N_ComputeBasisFunctions" :lisp "srf-nurb--st-compute-basis-functions" :file "NL_SrfNurb.cpp" :line 51 :lines 22 :wave 0 :status :done :static t
   :callees ())
 (:c "N_ComputeDerivatives" :lisp "srf-nurb--st-compute-derivatives" :file "NL_SrfNurb.cpp" :line 190 :lines 11 :wave 0 :status :done :static t
   :callees ())
 (:c "N_SrfBasisDerivs" :lisp "srf-basis-derivs" :file "NL_SrfNurb.cpp" :line 259 :lines 157 :wave 0 :status :done
   :callees ("N_BasisFindSpan" "N_ComputeBasisFunctions" "N_ComputeDerivatives" "N_ErrSet" "N_FloatOpIsBad" "N_IsSrfRat"))
 (:c "N_SrfEvalPt" :lisp "srf-eval-pt" :file "NL_SrfNurb.cpp" :line 1533 :lines 27 :wave 0 :status :done :root t
   :callees ())
 (:c "N_SrfEvalPtCrvOnSrf" :lisp "srf-eval-pt-crv-on-srf" :file "NL_SrfNurb.cpp" :line 1602 :lines 19 :wave 0 :status :port
   :callees ("N_ClampSrfAtParams" "N_CrvEval" "N_EndNurbs" "N_InitNurbs" "N_PtToXYZ" "N_SrfEvalPt"))
 (:c "N_ComputeDerivativeMatrixSpec" :lisp "compute-derivative-matrix-spec" :file "NL_SrfNurb.cpp" :line 1699 :lines 44 :wave 0 :status :done
   :callees ())
 (:c "N_ComputeDerivativeMatrix" :lisp "compute-derivative-matrix" :file "NL_SrfNurb.cpp" :line 1786 :lines 19 :wave 0 :status :done
   :callees ("N_ComputeDerivativeMatrixSpec"))
 (:c "N_SrfDerivs" :lisp "srf-derivs" :file "NL_SrfNurb.cpp" :line 1901 :lines 36 :wave 0 :status :done
   :callees ("N_ComputeDerivativeMatrix" "N_SrfBasisDerivs"))
 (:c "N_SrfGetBBox" :lisp "srf-get-b-box" :file "NL_SrfNurb.cpp" :line 1997 :lines 16 :wave 0 :status :done :root t
   :callees ("N_ENetGetBBox" "N_EndNurbs" "N_InitNurbs" "N_SrfGetArraySizes" "N_SrfGetENet"))
 (:c "N_SrfMinMaxWeightPosVectors" :lisp "srf-min-max-weight-pos-vectors" :file "NL_SrfNurb.cpp" :line 2188 :lines 48 :wave 0 :status :done
   :callees ("N_CPtGetW" "N_CPtToPtEuclid" "N_EndNurbs" "N_InitNurbs" "N_IsSrfRat" "N_SrfGetCPts" "N_VectorMagnitude"))
 (:c "N_SrfEvalPtPtDerivNormal" :lisp "srf-eval-pt-pt-deriv-normal" :file "NL_SrfNurb.cpp" :line 2310 :lines 14 :wave 0 :status :done :root t
   :callees ("N_AllocPt2dArray" "N_EndNurbs" "N_InitNurbs" "N_SrfEvalPtPtDerivNormalFast"))
 (:c "N_SrfEvalPtPtDerivNormalFast" :lisp "srf-eval-pt-pt-deriv-normal-fast" :file "NL_SrfNurb.cpp" :line 2382 :lines 36 :wave 0 :status :done
   :callees ("N_CopyPt" "N_SrfDerivs" "N_SrfEvalPtPtDerivNormalPole" "N_VectorCopy" "N_VectorCross" "N_VectorMagnitude" "N_VectorNormalizeRef"))
 (:c "N_SrfTranslate" :lisp "srf-translate" :file "NL_SrfNurb.cpp" :line 2873 :lines 13 :wave 0 :status :done
   :callees ("N_SrfGetCPts" "N_TranslateCPt"))
 (:c "N_SrfType" :lisp "srf-type" :file "NL_SrfNurb.cpp" :line 3769 :lines 576 :wave 0 :status :done
   :callees ("N_AllocCrvAndArrays" "N_AllocPt1dArray" "N_AllocReal1dArray" "N_CPtToPtEuclid" "N_Combine4Pts" "N_CreateLineStartDirVector" "N_CrvCopy" "N_CrvGetType" "N_CrvInitArrays" "N_CrvIsDegen" "N_CrvIsPlanar" "N_CrvLineFromPtAndVector" "N_CrvSetSizeIndices" "N_CrvTranslate" "N_CrvsAreCoincident" "N_DistPtInfLine" "N_DistPtLineSeg" "N_DistPtPt" "N_DistSignedPtPlane" "N_EndNurbs" "N_InitNurbs" "N_ProjectPtLine" "N_PtToXYZ" "N_SrfExtractIsoCrv" "N_SrfGetCPts" "N_SrfGetDegrees" "N_SrfGetKnots" "N_VectorCopy" "N_VectorCross" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude" "N_VectorNormalizeRef" "N_VectorReverse" "N_VectorScale" "N_VectorSum"))
 (:c "N_SrfReparam" :lisp "srf-reparam" :file "NL_SrfNurb.cpp" :line 4657 :lines 10 :wave 0 :status :done
   :callees ("N_CreateRectangle" "N_SrfGetParameterBounds" "N_SrfReparamToInterval"))
 (:c "N_SrfEvalPtPtDerivNormalPole" :lisp "srf-eval-pt-pt-deriv-normal-pole" :file "NL_SrfNurb.cpp" :line 5304 :lines 14 :wave 0 :status :done
   :callees ("N_AllocPt2dArray" "N_EndNurbs" "N_InitNurbs" "N_SrfEvalPtNormalAtPole"))
 (:c "N_SrfEvalPtNormalAtPole" :lisp "srf-eval-pt-normal-at-pole" :file "NL_SrfNurb.cpp" :line 5376 :lines 59 :wave 0 :status :done
   :callees ("N_CopyPt" "N_SrfDerivs" "N_SrfGetParameterBounds" "N_VectorCopy" "N_VectorCross" "N_VectorMagnitude" "N_VectorNormalizeRef" "N_VectorScale"))
 (:c "N_SrfMaxDiagDistBBox" :lisp "srf-max-diag-dist-b-box" :file "NL_SrfNurb.cpp" :line 5516 :lines 14 :wave 0 :status :done
   :callees ("N_BBoxGetDiagonal" "N_EndNurbs" "N_InitNurbs" "N_SrfGetBBox"))
 (:c "N_SrfOffset" :lisp "srf-offset" :file "NL_SrfOffset.cpp" :line 678 :lines 409 :wave 7 :status :done :root t
   :callees ("N_AllocInt1dArray" "N_AllocInt2dArray" "N_AllocKnotVectorAndArray" "N_AllocPt2dArray" "N_AllocReal1dArray" "N_CPtToPtAndW" "N_CreateRectangle" "N_CreateRevolvedSrf" "N_CreateSrfCornerPts" "N_CreateSrfExtrudeCrv" "N_CrvEvalTangent" "N_CrvGetParamBounds" "N_CrvInitArrays" "N_CrvIsPlanar" "N_CrvOffsetPtSampling" "N_EndNurbs" "N_FitCalcSrfParamValues" "N_FitCrvCalcKnotVector" "N_FitSrfLstSqApprox" "N_FitSrfToPtsKnots" "N_FloatOpIsBad" "N_FreePt2dArray" "N_InitNurbs" "N_PtsAreEqual" "N_SrfCompress" "N_SrfCopy" "N_SrfDecomposeToBez" "N_SrfEvalPtNormalAtPole" "N_SrfEvalPtPtDerivNormalFast" "N_SrfGetArraySizes" "N_SrfGetAverageLen" "N_SrfGetCPts" "N_SrfGetDegrees" "N_SrfGetGridPtsNormals" "N_SrfGetKnots" "N_SrfGetParameterBounds" "N_SrfInitArrays" "N_SrfOffsetGetMaxDeriv" "N_SrfRemoveAllKnots" "N_SrfReparamToInterval" "N_SrfTranslate" "N_SrfType" "N_SwapUV" "N_VectorCross" "N_VectorDiff" "N_VectorDot" "N_VectorMagnitude" "N_VectorNormalizeRef" "N_VectorPtAlongVector" "N_VectorScale" "N_VectorScaleRef" "N_Weight"))
 (:c "N_SrfInsertKnot" :lisp "srf-insert-knot" :file "NL_SrfTool.cpp" :line 68 :lines 169 :wave 1 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocReal2dArray" "N_AllocSrfArrays" "N_BasisFindSpanAndMult" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_FreeSrf" "N_InitNurbs" "N_KnotVectorIsEndParam" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_SrfSizeArrays"))
 (:c "N_SrfSplit" :lisp "srf-split" :file "NL_SrfTool.cpp" :line 403 :lines 183 :wave 1 :status :done :root t
   :callees ("N_BasisFindSpanAndMult" "N_Combine2CPts" "N_CopyCPt" "N_KnotVectorIsEndParam" "N_SrfGetCPtsDegreesAndKnots" "N_SrfSizeArrays"))
 (:c "N_SrfDecomposeToBez" :lisp "srf-decompose-to-bez" :file "NL_SrfTool.cpp" :line 870 :lines 170 :wave 1 :status :done
   :callees ("N_Alloc2dArraySrfPtrsParameters" "N_AllocCPt2dArray" "N_AllocReal1dArray" "N_BasisGetSpanCount" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_InitNurbs" "N_SrfGetCPts" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors"))
 (:c "N_SrfExtractIsoCrv" :lisp "srf-extract-iso-crv" :file "NL_SrfTool.cpp" :line 1173 :lines 152 :wave 1 :status :done :root t
   :callees ("N_BasisFindSpanAndMult" "N_Combine2CPts" "N_CopyCPt" "N_CrvSizeArrays" "N_SrfGetCPtsDegreesAndKnots"))
 (:c "N_SrfInsertKnots" :lisp "srf-insert-knots" :file "NL_SrfTool.cpp" :line 1466 :lines 197 :wave 1 :status :done
   :callees ("N_AllocSrfArrays" "N_BasisFindSpan" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_FreeSrf" "N_InitNurbs" "N_KnotVectorGetKnots" "N_KnotVectorIsEndParam" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_SrfSizeArrays"))
 (:c "N_SrfExtractPatch" :lisp "srf-extract-patch" :file "NL_SrfTool.cpp" :line 1814 :lines 195 :wave 1 :status :done
   :callees ("N_AllocCPt2dArray" "N_AllocSrfArrays" "N_BasisFindSpanAndMult" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_FreeSrf" "N_InitNurbs" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_SrfReparamToInterval" "N_SrfSizeArrays"))
 (:c "N_SrfRemoveKnotConditional" :lisp "srf-remove-knot-conditional" :file "NL_SrfTool.cpp" :line 2191 :lines 409 :wave 1 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_BasisFindGlobalMax" "N_BasisFindSpanAndMult" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_DistCptCptHomo" "N_EndNurbs" "N_InitNurbs" "N_IsSrfRat" "N_KnotVectorIsEndParam" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_SrfMinMaxWeightPosVectors" "N_SrfSetSizeIndices" "N_SrfSizeArrays"))
 (:c "N_SrfRemoveOneKnot" :lisp "srf-remove-one-knot" :file "NL_SrfTool.cpp" :line 2859 :lines 134 :wave 1 :status :done
   :callees ("N_AllocCPt1dArray" "N_AllocReal1dArray" "N_Combine2CPts" "N_CopyCPt" "N_DistCptCptHomo" "N_EndNurbs" "N_InitNurbs" "N_SrfGetCPtsDegreesAndKnots"))
 (:c "N_SrfRemoveAllKnots" :lisp "srf-remove-all-knots" :file "NL_SrfTool.cpp" :line 3114 :lines 694 :wave 1 :status :done :root t
   :callees ("N_AllocCPt2dArray" "N_AllocInt1dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_BasisFindAllSpanMaxima" "N_CPtGetW" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_InitNurbs" "N_IsSrfRat" "N_SrfCompress" "N_SrfCopy" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetCPtsKnotVectorAndKnots" "N_SrfMinMaxWeightPosVectors" "N_SrfRemoveOneKnot" "N_SrfSetSizeIndices" "N_SrfSizeArrays"))
 (:c "N_SrfElevateDegree" :lisp "srf-elevate-degree" :file "NL_SrfTool.cpp" :line 4224 :lines 424 :wave 1 :status :done
   :callees ("N_AllocCPt2dArray" "N_AllocReal1dArray" "N_AllocReal2dArray" "N_AllocSrfArrays" "N_BasisGetSpanCount" "N_BezGetDegreeElevationMatrix" "N_BezSrfElevateDegree" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_FreeSrf" "N_InitNurbs" "N_InitRealMatrix" "N_KnotVectorIsValid" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfGetKnotVectors" "N_SrfSizeArrays"))
 (:c "N_MakeSrfsCompatibleUV" :lisp "make-srfs-compatible-uv" :file "NL_SrfTool.cpp" :line 6107 :lines 4 :wave 1 :status :done :root t
   :callees ("N_MakeSrfsCompatible"))
 (:c "N_MakeSrfsCompatible" :lisp "make-srfs-compatible" :file "NL_SrfTool.cpp" :line 6151 :lines 111 :wave 1 :status :done
   :callees ("N_Alloc1dArrayKnotVectPtrs" "N_EndNurbs" "N_GetCompatibleKnotArray" "N_InitNurbs" "N_KnotVectorGetKnots" "N_MakeSrfsRatCompatible" "N_SrfElevateDegree" "N_SrfGetDegrees" "N_SrfGetKnotVectors" "N_SrfGetParameterBounds" "N_SrfInsertKnots" "N_SrfReparam" "N_SrfReparamToInterval"))
 (:c "N_SrfJoin" :lisp "srf-join" :file "NL_SrfTool.cpp" :line 6504 :lines 186 :wave 1 :status :done :root t
   :callees ("N_AllocSrfArrays" "N_CPtGetW" "N_CPtSetW" "N_Combine2CPts" "N_CopyCPt" "N_EndNurbs" "N_FreeSrf" "N_InitNurbs" "N_IsSrfRat" "N_SrfElevateDegree" "N_SrfGetCPtsAndKnots" "N_SrfGetCPtsDegreesAndKnots" "N_SrfInitArrays" "N_SrfRemoveKnotConditional" "N_SrfSizeArrays"))
 (:c "N_TessTrimmedSrf" :lisp "tess-trimmed-srf" :file "NL_Tessellate.cpp" :line 248 :lines 729 :wave 2 :status :port :root t
   :callees ("N_AllocInt1dArray" "N_AllocIntPtr1dArray" "N_AllocReal1dArray" "N_AllocRealPtr1dArray" "N_AllocRealPtr3dArray" "N_CPtToPtEuclid" "N_CrvDecomposeContinuity" "N_CrvEval" "N_CrvGetCPts" "N_CrvGetParamBounds" "N_DistCPolygon" "N_EndNurbs" "N_FloatOpIsBad" "N_FreeReal1dArray" "N_InitNurbs" "N_PtToXYZ" "N_Realloc1dFlagArray" "N_Realloc1dRealArray" "N_SrfGetParameterBounds" "N_SrfIsClosed" "N_TessSrfArea" "ST_TessIsOnBoundary" "ST_TessPolygonizeTrimCrv" "ST_TessRegionOfRectangles" "ST_TessSubdivideSrf"))
 (:c "ST_TessIsOnBoundary" :lisp "tessellate--st-tess-is-on-boundary" :file "NL_Tessellate.cpp" :line 1196 :lines 55 :wave 2 :status :port :static t
   :callees ("N_CPtToPtEuclid" "N_CrvGetCPts" "N_CrvIsLine" "N_PtToXYZ"))
 (:c "N_TessGetTriangleEdges" :lisp "tess-get-triangle-edges" :file "NL_Tessellate.cpp" :line 1317 :lines 47 :wave 2 :status :port :root t
   :callees ("N_AllocInt1dArray" "N_EndNurbs" "N_InitNurbs"))
 (:c "N_TessGetTriangles" :lisp "tess-get-triangles" :file "NL_Tessellate.cpp" :line 1519 :lines 54 :wave 2 :status :port :root t
   :callees ("N_AllocInt1dArray" "N_EndNurbs" "N_InitNurbs"))
 (:c "ST_TessIsInsideLoop" :lisp "tessellate--st-tess-is-inside-loop" :file "NL_Tessellate.cpp" :line 1657 :lines 42 :wave 2 :status :port :static t
   :callees ("ST_TessIsInsidePolygon"))
 (:c "ST_TessIsInsidePolygon" :lisp "tessellate--st-tess-is-inside-polygon" :file "NL_Tessellate.cpp" :line 1761 :lines 229 :wave 2 :status :port :static t
   :callees ("N_LineSegs2dAreIntersecting"))
 (:c "ST_TessIsInsideTrim" :lisp "tessellate--st-tess-is-inside-trim" :file "NL_Tessellate.cpp" :line 2106 :lines 46 :wave 2 :status :port :static t
   :callees ("ST_TessIsInsidePolygon"))
 (:c "ST_TessFindBBoxIndices" :lisp "tessellate--st-tess-find-b-box-indices" :file "NL_Tessellate.cpp" :line 2375 :lines 25 :wave 2 :status :port :static t
   :callees ("N_EndNurbs" "N_InitNurbs"))
 (:c "ST_TessRegionOfRectangles" :lisp "tessellate--st-tess-region-of-rectangles" :file "NL_Tessellate.cpp" :line 2984 :lines 734 :wave 2 :status :port :static t
   :callees ("N_AllocFlag1dArray" "N_AllocInt1dArray" "N_AllocInt2dArray" "N_AllocIntPtr1dArray" "N_AllocIntPtr2dArray" "N_EndNurbs" "N_FreeInt1dArray" "N_FreeIntPtr1dArray" "N_InitNurbs" "N_Realloc1dFlagArray" "N_Realloc1dIntArray" "N_Realloc2dIntArray" "N_Realloc2dIntPtrArray" "ST_TessBBoxXTrim" "ST_TessCosine" "ST_TessEdgeTouch" "ST_TessFindBBoxIndices" "ST_TessInitArray" "ST_TessIsEdgeOnBoundary" "ST_TessIsInsideLoop" "ST_TessLineXPolygon" "ST_TessPtInSpan" "ST_TessSideTest"))
 (:c "ST_TessInitArray" :lisp "tessellate--st-tess-init-array" :file "NL_Tessellate.cpp" :line 4048 :lines 6 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessSideTest" :lisp "tessellate--st-tess-side-test" :file "NL_Tessellate.cpp" :line 4060 :lines 9 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessPtInSpan" :lisp "tessellate--st-tess-pt-in-span" :file "NL_Tessellate.cpp" :line 4076 :lines 25 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessCosine" :lisp "tessellate--st-tess-cosine" :file "NL_Tessellate.cpp" :line 4113 :lines 19 :wave 2 :status :port :static t
   :callees ("N_FloatOpIsBad"))
 (:c "ST_TessIsEdgeOnBoundary" :lisp "tessellate--st-tess-is-edge-on-boundary" :file "NL_Tessellate.cpp" :line 4146 :lines 22 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessEdgeTouch" :lisp "tessellate--st-tess-edge-touch" :file "NL_Tessellate.cpp" :line 4177 :lines 55 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessLineXPolygon" :lisp "tessellate--st-tess-line-x-polygon" :file "NL_Tessellate.cpp" :line 4256 :lines 72 :wave 2 :status :port :static t
   :callees ())
 (:c "N_TessSrfArea" :lisp "tess-srf-area" :file "NL_Tessellate.cpp" :line 4389 :lines 50 :wave 2 :status :port
   :callees ("N_AllocPt1dArray" "N_CPtToPtEuclid" "N_DistPtPt" "N_EndNurbs" "N_InitNurbs" "N_SrfGetCPts"))
 (:c "ST_IsectLineWith2dPolygon" :lisp "tessellate--st-isect-line-with-2d-polygon" :file "NL_Tessellate.cpp" :line 4538 :lines 551 :wave 2 :status :port :static t
   :callees ("N_AllocFlag1dArray" "N_AllocInt1dArray" "N_AllocReal1dArray" "N_EndNurbs" "N_InitNurbs" "ST_ReallocArrays" "ST_TessLineXLine2d" "ST_TessOverlap" "ST_TessSideTest2d"))
 (:c "ST_TessSideTest2d" :lisp "tessellate--st-tess-side-test-2d" :file "NL_Tessellate.cpp" :line 5280 :lines 13 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessOverlap" :lisp "tessellate--st-tess-overlap" :file "NL_Tessellate.cpp" :line 5301 :lines 47 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_ReallocArrays" :lisp "tessellate--st-realloc-arrays" :file "NL_Tessellate.cpp" :line 5360 :lines 34 :wave 2 :status :port :static t
   :callees ("N_Realloc1dFlagArray" "N_Realloc1dIntArray" "N_Realloc1dRealArray"))
 (:c "ST_TessLineXLine2d" :lisp "tessellate--st-tess-line-x-line-2d" :file "NL_Tessellate.cpp" :line 5454 :lines 52 :wave 2 :status :port :static t
   :callees ("N_FloatOpIsBad"))
 (:c "ST_TessCheckTrim" :lisp "tessellate--st-tess-check-trim" :file "NL_Tessellate.cpp" :line 5570 :lines 60 :wave 2 :status :port :static t
   :callees ("N_BezInterp4Pts" "N_CPtToPtEuclid" "N_ClampSrfAtParams" "N_Combine4Pts" "N_CrvEval" "N_DistPtLineSeg" "N_DistPtPt" "N_EndNurbs" "N_InitNurbs" "N_PtToXYZ" "N_SrfEvalPt"))
 (:c "ST_TessBBoxXTrim" :lisp "tessellate--st-tess-b-box-x-trim" :file "NL_Tessellate.cpp" :line 5751 :lines 1311 :wave 2 :status :port :static t
   :callees ("N_AllocFlag1dArray" "N_AllocFlagPtr1dArray" "N_AllocInt1dArray" "N_AllocIntPtr1dArray" "N_AllocReal1dArray" "N_AllocRealPtr1dArray" "N_BBox2dCalcOverlap" "N_EndNurbs" "N_InitNurbs" "ST_IsPtInPoly" "ST_IsPtInPolygon" "ST_IsectLineWith2dPolygon" "ST_TessAppendPolygonVertex" "ST_TessAppendPt" "ST_TessBuffer" "ST_TessFindCrossingPt" "ST_TessFindUnused" "ST_TessIsPtBBoxCorner" "ST_TessIsPtInBBox"))
 (:c "ST_TessAppendPt" :lisp "tessellate--st-tess-append-pt" :file "NL_Tessellate.cpp" :line 7568 :lines 29 :wave 2 :status :port :static t
   :callees ("N_Realloc1dRealArray"))
 (:c "ST_TessBuffer" :lisp "tessellate--st-tess-buffer" :file "NL_Tessellate.cpp" :line 7615 :lines 47 :wave 2 :status :port :static t
   :callees ("N_Realloc1dIntArray"))
 (:c "ST_TessFindUnused" :lisp "tessellate--st-tess-find-unused" :file "NL_Tessellate.cpp" :line 7685 :lines 137 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessAppendPolygonVertex" :lisp "tessellate--st-tess-append-polygon-vertex" :file "NL_Tessellate.cpp" :line 7841 :lines 22 :wave 2 :status :port :static t
   :callees ("N_Realloc1dIntArray"))
 (:c "ST_TessIsPtInBBox" :lisp "tessellate--st-tess-is-pt-in-b-box" :file "NL_Tessellate.cpp" :line 7878 :lines 8 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessIsPtBBoxCorner" :lisp "tessellate--st-tess-is-pt-b-box-corner" :file "NL_Tessellate.cpp" :line 7893 :lines 23 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_IsPtInPolygon" :lisp "tessellate--st-is-pt-in-polygon" :file "NL_Tessellate.cpp" :line 7925 :lines 24 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_IsPtInPoly" :lisp "tessellate--st-is-pt-in-poly" :file "NL_Tessellate.cpp" :line 7960 :lines 25 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessFindCrossingPt" :lisp "tessellate--st-tess-find-crossing-pt" :file "NL_Tessellate.cpp" :line 7996 :lines 179 :wave 2 :status :port :static t
   :callees ())
 (:c "ST_TessPolygonizeTrimCrv" :lisp "tessellate--st-tess-polygonize-trim-crv" :file "NL_Tessellate.cpp" :line 8265 :lines 277 :wave 2 :status :port :static t
   :callees ("N_AllocReal1dArray" "N_BBoxGetDiagonal" "N_CPtToPtEuclid" "N_ClampSrfAtParams" "N_CopyPt" "N_CrvDecomposeContinuity" "N_CrvEval" "N_CrvGetBBox" "N_CrvGetCPts" "N_CrvGetDegree" "N_CrvGetKnots" "N_DistPtPt" "N_EndNurbs" "N_FloatOpIsBad" "N_InitNurbs" "N_PtToXYZ" "N_Realloc1dRealArray" "N_SrfEvalPtCrvOnSrf" "N_SrfGetCPts" "N_SrfGetParameterBounds" "N_VectorDiffCPts" "N_VectorDot" "ST_TessCheckTrim"))
 (:c "ST_TessSubdivideSrf" :lisp "tessellate--st-tess-subdivide-srf" :file "NL_Tessellate.cpp" :line 8749 :lines 441 :wave 2 :status :port :static t
   :callees ("N_AllocFlag1dArray" "N_AllocFlag2dArray" "N_AllocInt1dArray" "N_AllocInt2dArray" "N_AllocReal1dArray" "N_AllocSrfArrays" "N_CreateRectangle" "N_EndNurbs" "N_InitNurbs" "N_InitRealHash" "N_Realloc1dFlagArray" "N_Realloc1dIntArray" "N_Realloc1dRealArray" "N_RetrieveRealHash" "N_SearchInsertReal" "N_SrfExtractPatch" "N_SrfGetArraySizes" "N_SrfGetDegrees" "N_SrfIsFlat" "N_SrfReparamToInterval" "N_SrfSetSizeIndices" "N_UpdateRealHash" "ST_TessIsInsideTrim"))
 (:c "N_InitNurbs" :lisp "init-nurbs" :file "NL_Util.cpp" :line 52 :lines 55 :wave 0 :status :vanish :root t
   :callees ("N_FileClose" "N_FileOpen"))
 (:c "N_EndNurbs" :lisp "end-nurbs" :file "NL_Util.cpp" :line 134 :lines 453 :wave 0 :status :vanish :root t
   :callees ("N_FileClose" "N_FileOpen" "N_Free"))
 (:c "N_Malloc" :lisp "malloc" :file "NL_Util.cpp" :line 791 :lines 15 :wave 0 :status :vanish
   :callees ("N_ErrSet"))
 (:c "N_Free" :lisp "free" :file "NL_Util.cpp" :line 936 :lines 8 :wave 0 :status :vanish
   :callees ())
 (:c "N_AllocFlag1dArray" :lisp "alloc-flag-1d-array" :file "NL_Util.cpp" :line 1013 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocFlag2dArray" :lisp "alloc-flag-2d-array" :file "NL_Util.cpp" :line 1083 :lines 57 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocFlagPtr1dArray" :lisp "alloc-flag-ptr-1d-array" :file "NL_Util.cpp" :line 1197 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocInt1dArray" :lisp "alloc-int-1d-array" :file "NL_Util.cpp" :line 1267 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocInt2dArray" :lisp "alloc-int-2d-array" :file "NL_Util.cpp" :line 1337 :lines 57 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocIntPtr1dArray" :lisp "alloc-int-ptr-1d-array" :file "NL_Util.cpp" :line 1451 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocIntPtr2dArray" :lisp "alloc-int-ptr-2d-array" :file "NL_Util.cpp" :line 1521 :lines 57 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocReal1dArray" :lisp "alloc-real-1d-array" :file "NL_Util.cpp" :line 1635 :lines 24 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocReal2dArray" :lisp "alloc-real-2d-array" :file "NL_Util.cpp" :line 1705 :lines 57 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocReal3dArray" :lisp "alloc-real-3d-array" :file "NL_Util.cpp" :line 1819 :lines 96 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocRealPtr1dArray" :lisp "alloc-real-ptr-1d-array" :file "NL_Util.cpp" :line 1979 :lines 24 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocRealPtr3dArray" :lisp "alloc-real-ptr-3d-array" :file "NL_Util.cpp" :line 2163 :lines 24 :wave 0 :status :vanish
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocPt1dArray" :lisp "alloc-pt-1d-array" :file "NL_Util.cpp" :line 2447 :lines 24 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocPt2dArray" :lisp "alloc-pt-2d-array" :file "NL_Util.cpp" :line 2518 :lines 57 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocCPt1dArray" :lisp "alloc-c-pt-1d-array" :file "NL_Util.cpp" :line 2793 :lines 24 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_AllocCPt2dArray" :lisp "alloc-c-pt-2d-array" :file "NL_Util.cpp" :line 2866 :lines 57 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_XYZTo2dCPtArray" :lisp "xyz-to-2d-c-pt-array" :file "NL_Util.cpp" :line 3327 :lines 16 :wave 0 :status :done
   :callees ("N_AllocCPt2dArray" "N_CPtFromWxWyWz"))
 (:c "N_AllocPtPtr2dArray" :lisp "alloc-pt-ptr-2d-array" :file "NL_Util.cpp" :line 4029 :lines 57 :wave 0 :status :vanish :root t
   :callees ("N_ErrSet" "N_Free" "N_Malloc"))
 (:c "N_FileOpen" :lisp "file-open" :file "NL_Util.cpp" :line 4276 :lines 13 :wave 0 :status :vanish
   :callees ("N_ErrSet"))
 (:c "N_FileClose" :lisp "file-close" :file "NL_Util.cpp" :line 4375 :lines 4 :wave 0 :status :done
   :callees ())
 (:c "N_FreeFlag1dArray" :lisp "free-flag-1d-array" :file "NL_Util.cpp" :line 4910 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_FreeInt1dArray" :lisp "free-int-1d-array" :file "NL_Util.cpp" :line 4983 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_FreeInt2dArray" :lisp "free-int-2d-array" :file "NL_Util.cpp" :line 5055 :lines 40 :wave 0 :status :vanish
   :callees ("N_Free" "N_FreeInt1dArray"))
 (:c "N_FreeIntPtr1dArray" :lisp "free-int-ptr-1d-array" :file "NL_Util.cpp" :line 5132 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_FreeIntPtr2dArray" :lisp "free-int-ptr-2d-array" :file "NL_Util.cpp" :line 5205 :lines 40 :wave 0 :status :vanish
   :callees ("N_Free" "N_FreeIntPtr1dArray"))
 (:c "N_FreeReal1dArray" :lisp "free-real-1d-array" :file "NL_Util.cpp" :line 5282 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_FreeReal2dArray" :lisp "free-real-2d-array" :file "NL_Util.cpp" :line 5355 :lines 40 :wave 0 :status :vanish
   :callees ("N_Free" "N_FreeReal1dArray"))
 (:c "N_FreePt1dArray" :lisp "free-pt-1d-array" :file "NL_Util.cpp" :line 5509 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_FreePt2dArray" :lisp "free-pt-2d-array" :file "NL_Util.cpp" :line 5581 :lines 40 :wave 0 :status :vanish
   :callees ("N_Free" "N_FreePt1dArray"))
 (:c "N_FreeCPt1dArray" :lisp "free-c-pt-1d-array" :file "NL_Util.cpp" :line 5659 :lines 37 :wave 0 :status :vanish
   :callees ("N_Free"))
 (:c "N_FreeCPt2dArray" :lisp "free-c-pt-2d-array" :file "NL_Util.cpp" :line 5733 :lines 40 :wave 0 :status :vanish
   :callees ("N_Free" "N_FreeCPt1dArray"))
 (:c "N_Realloc1dFlagArray" :lisp "realloc-1d-flag-array" :file "NL_Util.cpp" :line 5812 :lines 17 :wave 0 :status :vanish
   :callees ("N_AllocFlag1dArray" "N_FreeFlag1dArray"))
 (:c "N_Realloc1dIntArray" :lisp "realloc-1d-int-array" :file "NL_Util.cpp" :line 5949 :lines 17 :wave 0 :status :vanish
   :callees ("N_AllocInt1dArray" "N_FreeInt1dArray"))
 (:c "N_Realloc2dIntArray" :lisp "realloc-2d-int-array" :file "NL_Util.cpp" :line 6083 :lines 21 :wave 0 :status :vanish
   :callees ("N_AllocInt2dArray" "N_FreeInt2dArray"))
 (:c "N_Realloc2dIntPtrArray" :lisp "realloc-2d-int-ptr-array" :file "NL_Util.cpp" :line 6154 :lines 21 :wave 0 :status :vanish
   :callees ("N_AllocIntPtr2dArray" "N_FreeIntPtr2dArray"))
 (:c "N_Realloc1dRealArray" :lisp "realloc-1d-real-array" :file "NL_Util.cpp" :line 6224 :lines 17 :wave 0 :status :vanish
   :callees ("N_AllocReal1dArray" "N_FreeReal1dArray"))
 (:c "N_Realloc2dRealArray" :lisp "realloc-2d-real-array" :file "NL_Util.cpp" :line 6291 :lines 21 :wave 0 :status :vanish
   :callees ("N_AllocReal2dArray" "N_FreeReal2dArray"))
 (:c "N_Realloc1dPtArray" :lisp "realloc-1d-pt-array" :file "NL_Util.cpp" :line 6361 :lines 17 :wave 0 :status :vanish
   :callees ("N_AllocPt1dArray" "N_CopyPt" "N_FreePt1dArray"))
 (:c "N_Realloc2dCPtArray" :lisp "realloc-2d-c-pt-array" :file "NL_Util.cpp" :line 6566 :lines 21 :wave 0 :status :vanish
   :callees ("N_AllocCPt2dArray" "N_CopyCPt" "N_FreeCPt2dArray"))
 (:c "N_PrintStack" :lisp "print-stack" :file "NL_Util.cpp" :line 6631 :lines 354 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_ErrSet" :lisp "err-set" :file "NL_Util.cpp" :line 8575 :lines 5 :wave 0 :status :vanish
   :callees ())
 (:c "N_ErrClear" :lisp "err-clear" :file "NL_Util.cpp" :line 8603 :lines 5 :wave 0 :status :done
   :callees ())
 (:c "N_CheckRealMatStorage" :lisp "check-real-mat-storage" :file "NL_Util.cpp" :line 8698 :lines 12 :wave 0 :status :done
   :callees ("N_ErrSet" "N_GetMaxIndexRealMatrix"))
 (:c "N_ErrGetType" :lisp "err-get-type" :file "NL_Util.cpp" :line 8857 :lines 4 :wave 0 :status :vanish :root t
   :callees ())
 (:c "N_CalcCircArcDerivs" :lisp "calc-circ-arc-derivs" :file "NL_Geom.cpp" :line 9131 :lines 32 :wave 0 :status :done
   :callees ("N_Combine2Pts" "N_TranslateSum2Pts"))
 (:c "N_ScaleCPtWithPtAndVector" :lisp "scale-c-pt-with-pt-and-vector" :file "NL_Math.cpp" :line 1160 :lines 27 :wave 0 :status :done
   :callees ())
 (:c "N_SrfRatToNonRat" :lisp "srf-rat-to-non-rat" :file "NL_SrfBasic.cpp" :line 4920 :lines 13 :wave 0 :status :done
   :callees ("N_CPtSetW" "N_SrfGetCPts"))
)
